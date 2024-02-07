package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.logging.Level;

import cern.colt.Arrays;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.subreads.NoopToSubreadsConverter;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.processing.DetectSnps;
import gedi.grand3.reads.ClippingData;
import gedi.grand3.reads.MismatchPerPositionStatistics;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.TargetCollection;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.EI;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;

public class Grand3SnpsAndClip<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3SnpsAndClip(Grand3ParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.nosnp);
		addInput(params.snp);
		addInput(params.snpConv);
		addInput(params.snpPval);
		addInput(params.strandness);
		addInput(params.experimentalDesignFile);
		addInput(params.clip5p1);
		addInput(params.clip5p2);
		addInput(params.clip3p1);
		addInput(params.clip3p2);
		addInput(params.noplot);
		addInput(params.targetCollection);
		addInput(params.blacklistSnp);
		
		addInput(params.prefix);
		addInput(params.debug);
		
		addOutput(params.snpFile);
		addOutput(params.strandnessFile);
		addOutput(params.mismatchPositionFile);
		addOutput(params.clipFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		int nthreads = getIntParameter(pind++);
		Genomic genomic = getParameter(pind++);
		GenomicRegionStorage<A> reads = getParameter(pind++);
		boolean nosnps = getBooleanParameter(pind++);
		String snps = getParameter(pind++);
		double conv = getDoubleParameter(pind++);
		double pvalCutoff = getDoubleParameter(pind++);
		Strandness strandness = getParameter(pind++);
		File expDesign = getParameter(pind++);
		int clip5p1 = getIntParameter(pind++);
		int clip5p2 = getIntParameter(pind++);
		int clip3p1 = getIntParameter(pind++);
		int clip3p2 = getIntParameter(pind++);
		boolean noplot = getBooleanParameter(pind++);
		TargetCollection targets = getParameter(pind++);
		boolean blacklistSnp = getBooleanParameter(pind++);
		
		String prefix = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		
		
		ExperimentalDesign design = ExperimentalDesign.fromTable(expDesign);
		
		context.getLog().info("Finding SNPs and collecting position statistics...");
		
		
		// Read count modes do not matter for snps (all are used, which calls the largest amount of snps)
		// for the posstat, only uniques are counted (thus, ints!) 
		targets = targets.create(ReadCountMode.Unique,ReadCountMode.Unique);
		targets.checkValid();
		
		context.getLog().info("Using the following categories for estimating clipping parameters: "+EI.wrap(targets.getCategories(c->c.useToInferClipping())).concat(","));
	
		ReadSource<A> source = new ReadSource<>(reads, null, Strandness.Unspecific, debug);
	
		context.getLog().info("Using "+StringUtils.removeFooter(source.getConverter().getClass().getSimpleName(),"ToSubreadsConverter")+" mode!");

		MismatchPerPositionStatistics posstat = source.getConverter().getClass().equals(NoopToSubreadsConverter.class)?null:new MismatchPerPositionStatistics(genomic,targets,design.getNumLibraries());
		if (posstat!=null) 
			posstat.setDebug(debug);

		boolean[] blacklist = new boolean[design.getCount()];
		if (blacklistSnp) {
			if (design.getTypes().length!=1) throw new RuntimeException("Blacklisting SNPs is only allowed when having a single label type!");
			
			IntArrayList blsamp = new IntArrayList();
			for (int c : design.getSamplesNotHaving(design.getTypes()[0])) {
				blsamp.add(c);
				for (int idx : design.getIndicesForSampleId(c))
					blacklist[idx]=true;
			}
			
			context.getLog().info("Blacklisting mismatches from samples "+blsamp.iterator().map(i->design.getSampleName(i)).concat(", "));

		}
		
		DetectSnps<A> algo = new DetectSnps<>(source, conv, pvalCutoff, blacklist);
		algo.setNthreads(nthreads);
		algo.process(context::getProgress, 
				targets,
				snps==null && !nosnps?getOutputFile(0).getPath():null, 
				posstat);
		
		if (nosnps) {
			EI.singleton("#").print(getOutputFile(0).getPath());
			context.getLog().info("Using no SNPs mode!");
		} else if(snps!=null) {
			EI.lines(snps).print(getOutputFile(0).getPath());
			context.getLog().info("SNPs overriden from "+snps);
		}

	
		context.getLog().info("Auto-Detecting sequencing mode: Sense:"+algo.getSense()+" Antisense:"+algo.getAntisense());
		Strandness inferredStrandness;
		if (algo.getSense()/2>algo.getAntisense()) {
			context.getLog().info("Detected strand-specific sequencing (Sense)");
			inferredStrandness = Strandness.Sense;
		} else if (algo.getAntisense()/2>algo.getSense()) {
			context.getLog().info("Detected strand-specific sequencing (Antisense)");
			inferredStrandness = Strandness.Antisense;
		} else {
			context.getLog().info("Detected strand-unspecific sequencing");
			inferredStrandness = Strandness.Unspecific;
		}
		if (strandness.equals(Strandness.AutoDetect)) strandness = inferredStrandness;
		else context.getLog().info("Overriden by command line: "+strandness.name());
		
		LineWriter s = getOutputWriter(1);
		s.writef("Sense\t%d\nAntisense\t%d\n", algo.getSense(),algo.getAntisense());
		s.writef("Inferred\t%s\nParameter\t%s\nUsed\t%s\n", inferredStrandness.toString(),getParameter(5).toString(), strandness.toString());
		s.close();
		
		int[] iclip1 = new int[2];
		int[] iclip2 = new int[2];

		if (posstat!=null) {
			context.getLog().info("Writing mismatches per position statistics...");
			posstat.write(targets.getCategories(),getOutputFile(2),strandness,reads.getMetaDataConditions());
		
			context.getLog().info("Inferring clipping parameters...");
			double[][] firstReadMatrix = null;
			double[][] secondReadMatrix = null;
			
			CompatibilityCategory[] cats = targets.getCategories(c->c.useToInferClipping());
			for (MetabolicLabelType l : design.getTypes()) {
				int[] samples = design.getLibraryIdsWithSamplesHaving(l);
				if (!strandness.equals(Strandness.Antisense)) {
					firstReadMatrix = ArrayUtils.cbind(firstReadMatrix,
							posstat.getMatrix(l.getGenomic(),l.getRead(),true,cats,true,samples));
					secondReadMatrix = ArrayUtils.cbind(secondReadMatrix,
							posstat.getMatrix(l.getGenomicReverse(),l.getReadReverse(),false,cats,false,samples));
				}
				
				if (!strandness.equals(Strandness.Sense)) {
					firstReadMatrix = ArrayUtils.cbind(firstReadMatrix,
							posstat.getMatrix(l.getGenomicReverse(),l.getReadReverse(),false,cats,true,samples));
					secondReadMatrix = ArrayUtils.cbind(secondReadMatrix,
							posstat.getMatrix(l.getGenomic(),l.getRead(),true,cats,false,samples));
				}
				
			}
		
//		RDataWriter rw = new RDataWriter(new FileOutputStream("trim.Rdata"));
//		rw.writeHeader();
//		rw.write("first", firstReadMatrix);
//		rw.write("second", secondReadMatrix);
//		rw.finish();

		
			iclip1 = ClippingData.inferClipping(firstReadMatrix,context.getLog());
			iclip2 = ClippingData.inferClipping(secondReadMatrix,context.getLog());
			
			if (secondReadMatrix==null || secondReadMatrix.length==0)
				context.getLog().info("Inferred: first read="+Arrays.toString(iclip1));
			else
				context.getLog().info("Inferred: first read="+Arrays.toString(iclip1)+" second read="+Arrays.toString(iclip2));
			if (clip5p1>=0 || clip5p2>=0 || clip3p1>=0 || clip3p2>=0)
				context.getLog().info("Overridden by parameter: first read=["+(clip5p1<0?"":clip5p1)+","+(clip3p1<0?"":clip3p1)+"] second read=["+(clip5p2<0?"":clip5p2)+","+(clip3p2<0?"":clip3p2)+"]");
			else if (clip5p1>posstat.getReadLength1()*0.2 || clip5p2>posstat.getReadLength2()*0.2 || clip3p1>posstat.getReadLength1()*0.2 || clip3p2>posstat.getReadLength2()*0.2) {
				context.getLog().severe("This is more than 20% of the read length. Check the position statistics, and specifiy clipping manually (-clip... x)!");
				System.exit(10);
			}

		}
		else {
			LineWriter dummy = getOutputWriter(2);
			dummy.writeLine("Mismatch statistics are not computed for Subread data!");
			dummy.close();
		}
		
		if (source.getConverter().getClass().equals(NoopToSubreadsConverter.class)) {
			context.getLog().info("Using no clipping by default for Noop...");
			clip5p1=clip3p1=clip5p2=clip3p2=0;
		}

		s = getOutputWriter(3);
		s.writef("Inferred 5p1\t%d\nInferred 5p2\t%d\nInferred 3p1\t%d\nInferred 3p2\t%d\n", iclip1[0],iclip2[0], iclip1[1],iclip2[1]);
		s.writef("Used 5p1\t%d\nUsed 5p2\t%d\nUsed 3p1\t%d\nUsed 3p2\t%d\n", 
				clip5p1<0?iclip1[0]:clip5p1,
				clip5p2<0?iclip2[0]:clip5p2,
				clip3p1<0?iclip1[1]:clip3p1,
				clip3p2<0?iclip2[1]:clip3p2
				);
		s.close();
		
		if (!noplot && posstat!=null) {
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".grand3.mismatchpos.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/grand3.mismatchpos.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		

		
		return null;
	}


	
}
