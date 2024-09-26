package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.Grand3Utils;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.SubreadCounterKNMatrices;
import gedi.grand3.processing.SubreadCountOverallMismatches;
import gedi.grand3.processing.SubreadCounter;
import gedi.grand3.processing.SubreadProcessor;
import gedi.grand3.reads.ClippingData;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.SnpData;
import gedi.grand3.targets.TargetCollection;
import gedi.util.functions.EI;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;

public class Grand3CollectSufficientStatistics<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3CollectSufficientStatistics(Grand3ParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.snpFile);
		addInput(params.strandnessFile);
		addInput(params.clipFile);
		addInput(params.experimentalDesignFile);
		addInput(params.noplot);
		addInput(params.targetCollection);
		
		addInput(params.prefix);
		addInput(params.debug);
		
		addOutput(params.knMatFile);
		addOutput(params.mismatchFile);
		addOutput(params.perrFile);
		addOutput(params.subreadSemanticFile);
		
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		int nthreads = getIntParameter(pind++);
		Genomic genomic = getParameter(pind++);
		GenomicRegionStorage<A> reads = getParameter(pind++);
		File snpFile = getParameter(pind++); 
		File strandnessFile = getParameter(pind++); 
		File clipFile = getParameter(pind++);
		File designFile = getParameter(pind++);
		boolean noplot = getBooleanParameter(pind++);
		TargetCollection targets = getParameter(pind++);
		
		String prefix = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		context.getLog().info("Collecting data for global parameter estimation");
		
		SnpData masked = new SnpData(snpFile);
		context.getLog().info("Read "+masked.size()+" SNPs!");
		Strandness strandness = Grand3Utils.getStrandness(strandnessFile);
		context.getLog().info("Strandness: "+strandness);
		ClippingData clipping = ClippingData.fromFile(clipFile);
		context.getLog().info("Clipping: "+clipping);
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);
		
//		if (strandness.equals(Strandness.Unspecific))
//			throw new RuntimeException("Strand unspecific libraries are not yet supported!");

		ReadSource<A> source = new ReadSource<>(reads, clipping, strandness, debug);
		Grand3Utils.writeSemantic(source.getConverter().getSemantic(),getOutputFile(3));
		
		targets = targets.create(ReadCountMode.Unique,ReadCountMode.Unique);
		targets.checkValid();
		context.getLog().info("Using the following categories for estimating global parameters: "+EI.wrap(targets.getCategories(c->c.useToEstimateGlobalParameters())).concat(","));
		
		SubreadCountOverallMismatches overallBarc = new SubreadCountOverallMismatches(targets.getNumCategories(),source.getConverter().getSemantic().length,design.getIndexToIndex());
		SubreadCountOverallMismatches overallSam = new SubreadCountOverallMismatches(targets.getNumCategories(),source.getConverter().getSemantic().length,design.getIndexToSampleId());
		SubreadCounterKNMatrices binom = new SubreadCounterKNMatrices(c->c.useToEstimateGlobalParameters(),design.getNumSamples(), design.getIndexToSampleId(),source.getConverter().getSemantic().length,design.getTypes());
		binom.setDebug(debug);
		
		ArrayList<SubreadCounter> counters = new ArrayList<>();
		counters.add(overallBarc);
		counters.add(overallSam);
		counters.add(binom);
		
//		EI.wrap(targets.getCategories())
//			.map(cat->new SubreadCounterKNMatrices(c->c.equals(cat), design.getNumSamples(), design.getIndexToSampleId(),source.getConverter().getSemantic().length,design.getTypes()))
//			.add(counters);
		
		
		// for the globalstat, only uniques are counted (thus, ints!) 
		SubreadProcessor<A> algo = new SubreadProcessor<>(genomic,source,masked,context.getLog());
		algo.setNthreads(nthreads);
		algo.setDebug(debug);
		algo.process(context::getProgress, 
				targets,counters.toArray(new SubreadCounter[0]));
//				overallBarc,overallSam,binom);
		
		
		context.getLog().info("Writing sufficient statistics");
		binom.write(getOutputFile(0), design);
		
//		for (int i=3; i<counters.size(); i++) {
//			SubreadCounterKNMatrices cbinom = (SubreadCounterKNMatrices) counters.get(i);
//			cbinom.write(new File(getOutputFile(0).getPath().replace("tsv.gz",targets.getCategories()[i-3]+".tsv.gz")), design);
//		}
		
		context.getLog().info("Writing quality metrics");
		overallBarc.write(getOutputFile(1), i->design.getFullName(i), i->design.getSampleName(i), targets);
		overallSam.writePerr(getOutputFile(2), i->design.getSampleNameForSampleIndex(i), targets);
		
		
		if (!noplot) {
			try {
				
				context.getLog().info("Running R scripts for plotting");
				
				RRunner r = new RRunner(prefix+".grand3.mismatchfreq.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/grand3.conversionfreq.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		return null;
	}


	
}
