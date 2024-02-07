package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.Grand3Utils;
import gedi.grand3.estimation.MismatchMatrix;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.TargetEstimationResult;
import gedi.grand3.estimation.TargetEstimationResult.ModelType;
import gedi.grand3.estimation.TargetEstimationResult.SingleEstimationResult;
import gedi.grand3.estimation.estimators.ModelEstimationMethod;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.processing.Resimulator;
import gedi.grand3.processing.Resimulator.ResimulatorModelType;
import gedi.grand3.processing.SubreadProcessor;
import gedi.grand3.reads.ClippingData;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.SnpData;
import gedi.grand3.targets.TargetCollection;
import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3Resimulate<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3Resimulate(Grand3ParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.targetCollection);
		addInput(params.snpFile);
		addInput(params.strandnessFile);
		addInput(params.clipFile);
		addInput(params.experimentalDesignFile);
		addInput(params.subreadSemanticFile);
		addInput(params.targetBinFile);
		addInput(params.modelFile);
		addInput(params.modelBinFile);
		addInput(params.estimMethod);
		addInput(params.resimModelType);
		addInput(params.resimSeed);
		addInput(params.perrFile);
		
		addInput(params.prefix);
		addInput(params.debug);
		
		setRunFlag(params.resim);
		
		addOutput(params.resimFile);
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
		File subreadFile = getParameter(pind++); 
		File targetFile = getParameter(pind++); 
		File modelFile = getParameter(pind++);
		File modelBinFile = getParameter(pind++);
		ArrayList<ModelEstimationMethod> estim = getParameters(pind++);
		ResimulatorModelType modelType = getParameter(pind++);
		long seed = getLongParameter(pind++);
		File mmFile = getParameter(pind++);
		
		String prefix = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		context.getLog().info("Processing targets");
		
		SnpData masked = new SnpData(snpFile);
		Strandness strandness = Grand3Utils.getStrandness(strandnessFile);
		ClippingData clipping = ClippingData.fromFile(clipFile);
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);
		String[] subreads = Grand3Utils.readSemantic(subreadFile);

		context.getLog().info("Read "+masked.size()+" SNPs!");
		context.getLog().info("Strandness: "+strandness);
		context.getLog().info("Clipping: "+clipping);
		
		MismatchMatrix mmMat = new MismatchMatrix(mmFile, design, subreads);
		
		ModelStructure[][][] models;
		if (modelFile.lastModified()>modelBinFile.lastModified()) {
			context.getLog().info("Reading model parameters from text file: "+estim.get(0));
			models = Grand3Utils.readModelsTsv(modelFile,design,subreads,estim.get(0));
		}
		else {
			context.getLog().info("Reading model parameters from binary file.");
			models = Grand3Utils.readModels(modelBinFile,design.getNumSamples());
		}

		context.getLog().info("Reading gene parameters "+getOutputFile(0));
		HashMap<String,Double>[][][] ntrs = readGeneParameters(targetFile,design,targets.getCategories(),modelType); 
		
		
		
		context.getLog().info("Writing resimulated (replace mismatches by: "+modelType+") subreads to CIT file "+getOutputFile(0));
		ReadSource<A> source = new ReadSource<>(reads, clipping, strandness, debug);
		
		Resimulator resim = ntrs==null?null:new Resimulator(genomic,modelType,models,ntrs,mmMat,design.getTypes(), seed);
		
		SubreadProcessor<A> algo = new SubreadProcessor<A>(genomic,source,masked);
		algo.setNthreads(nthreads);
		algo.setDebug(debug);
		algo.writeSubreads(context::getProgress, 
				targets,
				getOutputFile(0),
				design,
				subreads,
				resim);
		
		
		return null;
	}


	/**
	 * t,i,cat
	 * @param targetFile
	 * @param design
	 * @param subreads
	 * @return
	 * @throws IOException
	 */
	private HashMap<String, Double>[][][] readGeneParameters(File targetFile, ExperimentalDesign design, CompatibilityCategory[] categories, ResimulatorModelType model) throws IOException {
		if (model==ResimulatorModelType.None) return null;
		
		HashMap<String, CompatibilityCategory> name2cat = EI.wrap(categories).index(c->c.getName());
		
		PageFile pf = new PageFile(targetFile.getPath());
		
		int storedNumFeatures = pf.getInt();
		long storedNumCounts = pf.getLong();
		long[] storedNumNtr = new long[design.getTypes().length];
		for (int i=0; i<storedNumNtr.length; i++)
			storedNumNtr[i] = pf.getLong();

		HashMap<String, Double>[][][] re = new HashMap[design.getTypes().length][design.getCount()][categories[categories.length-1].id()+1];
		for (int t=0; t<re.length; t++)
			for (int i=0; i<design.getCount(); i++)
				for (int c=0; c<categories.length; c++)
					if (categories[c].useToEstimateTargetParameters())
						re[t][i][c] = new HashMap<String, Double>();
		
		
		for (TargetEstimationResult res : pf.ei(()->new TargetEstimationResult()).loop()) {
			if (res.isEmpty()) continue;
			
			CompatibilityCategory cat = name2cat.get(res.getCategory());
			if (cat.useToEstimateTargetParameters())
				for (int t=0; t<re.length; t++) {
					IntIterator it = res.iterate(t);
					while (it.hasNext()) {
						int i = it.nextInt();
						SingleEstimationResult val = res.get(t,i);
						
						re[t][i][cat.id()].put(res.getTarget().getData(), val.getModel(ModelType.valueOf(model.name())).getMap());
					}
				}
			
		}
		
		pf.close();
		
		return re;
	}

	
	
}
