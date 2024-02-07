package gedi.slam.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.regex.Pattern;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.slam.GeneData;
import gedi.slam.OptimNumericalIntegrationProportion;
import gedi.slam.SlamCollector;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.dataframe.DataFrame;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.tsv.formats.Csv;
import gedi.util.math.stat.factor.Factor;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class SlamTest extends GediProgram {

	
	
	public SlamTest(SlamParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.snpFile);
		addInput(params.locations);
		addInput(params.strandness);
		addInput(params.trim5p);
		addInput(params.trim3p);
		addInput(params.rateTable);
		addInput(params.no4sUpattern);
		
		setRunFlag(params.test);
		addOutput(params.testFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException {
		
		int nthreads = getIntParameter(0);
		Genomic genomic = getParameter(1);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(2);
		File snps = getParameter(3);
		GenomicRegionStorage<NameProvider> locations = getParameter(4);
		Strandness strandness = getParameter(5);
		int trim5p = getIntParameter(6);
		int trim3p = getIntParameter(7);
		File rateFile = getParameter(8);
		String pat = getParameter(9);
		
		
		
		Pattern no4sUPattern = Pattern.compile(pat,Pattern.CASE_INSENSITIVE);
		String[] conds = reads.getMetaDataConditions();
		int[] no4sUIndices = EI.along(conds).filterInt(i->no4sUPattern.matcher(conds[i]).find()).toIntArray();
		boolean[] no4sU = new boolean[conds.length];
		for (int i : no4sUIndices)
			no4sU[i] = true;
		
		context.getLog().info("Reading SNPs...");
		DataFrame snpDf = Csv.toDataFrame(snps.getAbsolutePath());
		
		MemoryIntervalTreeStorage<Void> masked = new MemoryIntervalTreeStorage<>(Void.class);
		if (snpDf.rows()>0)
			for (Factor f : (Factor[]) snpDf.getColumn(0).getRaw()) 
				masked.add(ImmutableReferenceGenomicRegion.parse(f.name()));

		context.getLog().info("Read "+snpDf.rows()+" SNPs!");
		
		
		context.getLog().info("Reading Rates...");
		ExtendedIterator<String[]> rit = EI.lines(rateFile.getPath()).map(s->StringUtils.split(s,'\t'));
		String[] cond=ArrayUtils.slice(rit.next(),1);
		
		double[] polda = new double[cond.length];  Arrays.fill(polda, -1);
		double[] pnewa = new double[cond.length];  Arrays.fill(pnewa, -1);
		double[] dpolda = new double[cond.length];  Arrays.fill(dpolda, -1);
		double[] dpnewa = new double[cond.length];  Arrays.fill(dpnewa, -1);
		
		rit.sideEffect(a->a[0].startsWith("single_old"),a->parseRates(polda,a))
			.sideEffect(a->a[0].startsWith("single_new"),a->parseRates(pnewa,a))
			.sideEffect(a->a[0].startsWith("double_old"),a->parseRates(dpolda,a))
			.sideEffect(a->a[0].startsWith("double_new"),a->parseRates(dpnewa,a))
			.drain();

		OptimNumericalIntegrationProportion[] vb = new OptimNumericalIntegrationProportion[cond.length];
		for (int i=0; i<vb.length; i++) {
			double pold = correct(polda[i]);
			double pnew = correct(pnewa[i]);
			double dpold = correct(dpolda[i]);
			double dpnew = correct(dpnewa[i]);
			
			System.out.println("For "+cond[i]+", using p_old = "+pold+"  p_new="+pnew+" dp_old="+dpold+" dp_new="+dpnew);
			
			vb[i] = new OptimNumericalIntegrationProportion(1, 1, pold,pnew,dpold,dpnew, 0.05, 0.95,true);
		}

		
		context.getLog().info("Testing reads (Strand: "+strandness.name()+") ...");
		
		
		
		SlamCollector collector = new SlamCollector(genomic, g->true, reads, masked, locations,strandness,trim5p,trim3p, ReadCountMode.Unique, ReadCountMode.Unique,no4sU,false,false,false,false);
//		collector.countNumis((g,r)->{
//			for (int i=0; i<reads.getMetaDataConditions().length; i++){
//				System.out.println(i+" "+reads.getMetaDataConditions()[i]+": "+r[i]);
//			}
//		});
//		collector.setVerbose(true);
		for (GeneData g : collector.collect("N").loop()) {
			System.out.println(g.getGene());
			for (int i=0; i<reads.getMetaDataConditions().length; i++){
				System.out.println(" "+reads.getMetaDataConditions()[i]+": "+g.getTotalConversions().getDouble(i)+"/"+g.getTotalCoverage().getDouble(i)+" "+g.getTotalDoubleHits().getDouble(i)+"/"+g.getTotalDoubleHitCoverage().getDouble(i)+"\t");
			}
		}
		
		
		return null;
	}

	
	private double correct(double d) {
		return Math.min(Math.max(d, 1E-6), 0.9);
	}


	private void parseRates(double[] rates, String[] a) {
		if (rates.length!=a.length-1) throw new RuntimeException("Rates file incompatible with reads!");
		for (int i=0; i<rates.length; i++)
			rates[i] = Double.parseDouble(a[i+1]);
	}
	
}
