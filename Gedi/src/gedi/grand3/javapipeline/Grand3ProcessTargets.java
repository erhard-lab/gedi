package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.logging.Level;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.Grand3Utils;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.TargetEstimationResult;
import gedi.grand3.estimation.TargetEstimator;
import gedi.grand3.estimation.estimators.ModelEstimationMethod;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.PseudobulkDefinition;
import gedi.grand3.knmatrix.SubreadCounterKNMatrices;
import gedi.grand3.knmatrix.SubreadCounterKNMatrixPerTarget;
import gedi.grand3.processing.SubreadProcessor;
import gedi.grand3.processing.SubreadStatistics;
import gedi.grand3.reads.ClippingData;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.SnpData;
import gedi.grand3.targets.TargetCollection;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.mapping.OneToManyMapping;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;

public class Grand3ProcessTargets<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3ProcessTargets(Grand3ParameterSet params, boolean hasMappedTarget,boolean hasTargetMixmat) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.targetCollection);
		addInput(params.snpFile);
		addInput(params.strandnessFile);
		addInput(params.clipFile);
		addInput(params.experimentalDesignFile);
		addInput(params.modelFile);
		addInput(params.modelBinFile);
		addInput(params.estimMethod);
		addInput(params.subreadSemanticFile);
		addInput(params.subreadsToUse);
		addInput(params.subreadCit);
		addInput(params.noplot);
		addInput(params.pseudobulkFile);
		addInput(params.pseudobulkName);
		addInput(params.pseudobulkMinimalPurity);
		addInput(params.targetMixmat);
		addInput(params.targetMergeTable);
		addInput(params.outputMixBeta);
		
		addInput(params.prefix);
		addInput(params.debug);
		
		if (hasTargetMixmat) {
			addOutput(params.targetMixmatFile);
		} else if (hasMappedTarget)
			addOutput(params.pseudobulkBinFile);
		else
			addOutput(params.targetBinFile);
		
		
		
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		int nthreads = getIntParameter(pind++);
		Genomic genomic = getParameter(pind++);
		GenomicRegionStorage<A> reads = getParameter(pind++);
		TargetCollection targets = getParameter(pind++);
		File snpFile = getParameter(pind++); 
		File strandnessFile = getParameter(pind++); 
		File clipFile = getParameter(pind++);
		File designFile = getParameter(pind++);
		File modelFile = getParameter(pind++);
		File modelBinFile = getParameter(pind++);
		ArrayList<ModelEstimationMethod> estim = getParameters(pind++);
		File subreadFile = getParameter(pind++);
		HashSet<String> useSubreads = new HashSet<>(getParameters(pind++));
		boolean sreadCit = getBooleanParameter(pind++);
		boolean noplot = getBooleanParameter(pind++);
		String pseudobulkFile = getParameter(pind++);
		String pseudobulkName = getParameter(pind++);
		double pseudobulkMinimalPurity = getDoubleParameter(pind++);
		String targetMixmat = getParameter(pind++);
		String targetMergeTab = getParameter(pind++);
		boolean writeMix = getBooleanParameter(pind++);
		
		String prefix = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		context.getLog().info("Processing targets");
		
		SnpData masked = new SnpData(snpFile);
		Strandness strandness = Grand3Utils.getStrandness(strandnessFile);
		ClippingData clipping = ClippingData.fromFile(clipFile);
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);
		String[] subreads = Grand3Utils.readSemantic(subreadFile);

		ModelStructure[][][] models;
		if (modelFile.lastModified()>modelBinFile.lastModified()) {
			context.getLog().info("Reading model parameters from text file: "+estim.get(0));
			models = Grand3Utils.readModelsTsv(modelFile,design,subreads,estim.get(0));
		}
		else {
			context.getLog().info("Reading model parameters from binary file.");
			models = Grand3Utils.readModels(modelBinFile,design.getNumSamples());
		}
		
		// which subreads to use
		boolean[] subreadsToUse = new boolean[subreads.length];
//		System.out.println(useSubreads.size());
//		for (String s : useSubreads)
//			System.out.println(s);
		if (useSubreads.size()==0)
			Arrays.fill(subreadsToUse, true);
		else
			for (int i=0; i<subreadsToUse.length; i++)
				subreadsToUse[i] = useSubreads.contains(i+"") || useSubreads.contains(subreads[i]);

		
		context.getLog().info("Read "+masked.size()+" SNPs!");
		context.getLog().info("Strandness: "+strandness);
		context.getLog().info("Clipping: "+clipping);
		context.getLog().info("Subreads used for parameter estimation: "+EI.seq(0,subreads.length).filterInt(i->subreadsToUse[i]).map(i->subreads[i]).concat(","));
		
		targets.checkValid();
		context.getLog().info("Using the following categories for estimating target parameters: "+EI.wrap(targets.getCategories(c->c.useToEstimateTargetParameters())).concat(","));
		
		int[][] targetMapping;
		int[] targetToSample;
		if (pseudobulkFile!=null) {
			context.getLog().info("Mapping output conditions using the file "+pseudobulkFile);
			PseudobulkDefinition psdef = new PseudobulkDefinition(pseudobulkFile,design,context.getLog(),pseudobulkMinimalPurity);
			targetMapping = psdef.getCellsToPseudobulk();
			targetToSample = psdef.getSampleForPseudobulks();
		}
		else {
			context.getLog().info("Using output conditions as in experimental design");
			targetMapping = new int[design.getCount()][];
			for (int i=0; i<targetMapping.length; i++)
				targetMapping[i] = new int[] {i};
			targetToSample = design.getIndexToSampleId();
		}
		if (EI.wrap(targetMapping).mapToInt(a->ArrayUtils.max(a)).max()+1!=targetToSample.length) throw new RuntimeException("Assertion failed!");
		context.getLog().info("Output conditions n="+targetToSample.length);
		
		ReadSource<A> source = new ReadSource<>(reads, clipping, strandness, debug);
		
		SubreadProcessor<A> algo = new SubreadProcessor<A>(genomic,source,masked,context.getLog());
		algo.setNthreads(nthreads);
		algo.setDebug(debug);
		
		if (targetMixmat!=null) {
			context.getLog().info("Will only compute parameters and mix matrix for "+targetMixmat+"...");
			SubreadCounterKNMatrices binom = new SubreadCounterKNMatrices(c->c.useToEstimateTargetParameters(),design.getNumSamples(), targetToSample,source.getConverter().getSemantic().length,design.getTypes());
			binom.setDebug(debug);
			algo.processTarget(context::getProgress,targetMixmat, targets, binom);
			binom.write(getOutputFile(0), design);
			
			return null;
		}
		
		context.getLog().info("Filling likelihood cache...");
		TargetEstimator targetEstimatorObject = new TargetEstimator(design, models, targetMapping, targetToSample, 0.01,context::getProgress,nthreads,writeMix);
		context.getLog().info("Done filling likelihood cache.");

		
		OneToManyMapping<String, String> mapper = new OneToManyMapping<String, String>();
		if (targetMergeTab!=null) {
			
			
			mapper = OneToManyMapping.fromFile(targetMergeTab, "merged", "name");
			context.getLog().info("Loaded merge table with "+mapper.getFromUniverse().size()+" entries.");
		}
		
		SubreadCounterKNMatrixPerTarget targetEstimator = new SubreadCounterKNMatrixPerTarget(targets.getCategories(),c->c.useToEstimateTargetParameters(), targetEstimatorObject, design.getTypes(),subreadsToUse,mapper);
		
		
		PageFileWriter bin = new PageFileWriter(getOutputFile(0).getPath());
		bin.putInt(0); // reserved
		bin.putLong(0); // reserved
		for (int i=0; i<design.getTypes().length; i++)
			bin.putLong(0); // reserved
		
		context.getLog().info("Estimating target parameters...");
		
		
		SubreadStatistics stat = new SubreadStatistics(source.getConverter().getSemantic().length, design.getIndexToSampleId());
		
		int num = 0;
		long nonzeroCounts = 0;
		long[] nonzeroNtrs = new long[design.getTypes().length];
		for (TargetEstimationResult res : algo.process(context::getProgress, 
				targets,mapper,
				targetEstimator,
				stat).loop()) 
			if (!res.isEmpty()) {
				res.serialize(bin);
				num++;
				nonzeroCounts+=res.getNumNonzeroCounts();
				for (int t=0; t<nonzeroNtrs.length; t++)
					nonzeroNtrs[t]+=res.getNumNonzeroNtrs(t);
			}
		
		bin.position(0);
		bin.putInt(num);
		bin.putLong(nonzeroCounts);
		for (int t=0; t<nonzeroNtrs.length; t++)
			bin.putLong(nonzeroNtrs[t]);
		
		bin.close();
		
		stat.writeLengths(new File(prefix+".reads.lengths.tsv"), i->design.getSampleNameForSampleIndex(i));
		stat.writeSubreads(new File(prefix+".reads.subreads.tsv"), i->design.getSampleNameForSampleIndex(i),source.getConverter().getSemantic());
		
		if (!noplot) {
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".grand3.readstat.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/grand3.readstat.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		return null;
	}

	
	
}
