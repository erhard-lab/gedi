package gedi.slam.javapipeline;


import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Pattern;

import cern.colt.bitvector.BitVector;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import jdistlib.Beta;
import jdistlib.Binomial;

public class SlamDetectSnps extends GediProgram {

	
	
	public SlamDetectSnps(SlamParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.snpConv);
		addInput(params.snpPval);
		addInput(params.strandness);
		addInput(params.no4sUpattern);
		addInput(params.newsnp);
		
		addInput(params.prefix);
		
		addOutput(params.snpFile);
		addOutput(params.strandnessFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException {
		
		int nthreads = getIntParameter(0);
		Genomic genomic = getParameter(1);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(2);
		double conv = getDoubleParameter(3);
		double pvalCutoff = getDoubleParameter(4);
		Strandness strandness = getParameter(5);
		String pat = getParameter(6);
		boolean newsnp = getBooleanParameter(7);
		
		Pattern no4sUPattern = Pattern.compile(pat,Pattern.CASE_INSENSITIVE);
		
		String[] cond = reads.getMetaDataConditions();
		int[] no4sUIndices = EI.along(cond).filterInt(i->no4sUPattern.matcher(cond[i]).find()).toIntArray();
		boolean[] no4sU = new boolean[cond.length];
		for (int i : no4sUIndices)
			no4sU[i] = true;
		
		
		if (newsnp)
			context.getLog().info("Finding SNPs (only using no4sU)...");
		else
			context.getLog().info("Finding SNPs...");
		
		AtomicInteger senseCounter = new AtomicInteger();
		AtomicInteger antisenseCounter = new AtomicInteger();
		genomic.getGenes().ei()
			.progress(context.getProgress(), (int)genomic.getGenes().size(), (r)->"In "+r.getData())
			.parallelized(nthreads, 5, ei->ei.map(l->newsnp?collectNew(pvalCutoff,reads,no4sU,l, senseCounter, antisenseCounter):collect(conv,pvalCutoff,reads,l, senseCounter, antisenseCounter)))
			.filter(s->s.length()>0)
			.print("Location\tCoverage\tMismatches\tP value",getOutputFile(0).getPath());
	
		context.getLog().info("Auto-Detecting sequencing mode: Sense:"+senseCounter.get()+" Antisense:"+antisenseCounter.get());
		if (senseCounter.get()>antisenseCounter.get()*2) {
			context.getLog().info("Detected strand-specific sequencing (Sense)");
			if (strandness.equals(Strandness.AutoDetect)) strandness = Strandness.Sense;
			else context.getLog().info("Overriden by command line: "+strandness.name());
		} else if (antisenseCounter.get()>senseCounter.get()*2) {
			context.getLog().info("Detected strand-specific sequencing (Antisense)");
			if (strandness.equals(Strandness.AutoDetect)) strandness = Strandness.Antisense;
			else context.getLog().info("Overriden by command line: "+strandness.name());
		} else {
			context.getLog().info("Detected strand-unspecific sequencing");
			if (strandness.equals(Strandness.AutoDetect)) strandness = Strandness.Unspecific;
			else context.getLog().info("Overriden by command line: "+strandness.name());
		}
		
		LineWriter sout = getOutputWriter(1);
		sout.writeLine(strandness.name());
		sout.close();
		
		return null;
	}


	public static String collect(double conv, double pvalCutoff, GenomicRegionStorage<AlignedReadsData> reads, ImmutableReferenceGenomicRegion<String> gene, AtomicInteger senseCounter, AtomicInteger antisenseCounter) {
		gene=gene.toMutable().transformRegion(r->r.extendBack(1000).extendFront(1000)).toImmutable();
		
		
		CoverageAlgorithm cov = new CoverageAlgorithm(gene).setExample(NumericArray.wrap(0.0));
		
		int[] co = new int[2];
		ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> rit = 
				reads.ei(gene).sideEffect(r->co[0]++)
				.chain(reads.ei(gene.toMutable().toOppositeStrand()).sideEffect(r->co[1]++));
		BitVector counted = new BitVector(300);
		
		HashMap<Integer,double[]> counter = new HashMap<>();
		for (ImmutableReferenceGenomicRegion<AlignedReadsData> r : rit.loop()) {
			if (r.getRegion().getTotalLength()>counted.size()) 
				counted = new BitVector(r.getRegion().getTotalLength());
			
			for (int d=0; d<r.getData().getDistinctSequences(); d++) {
				counted.clear();
				double c = r.getData().getTotalCountForDistinct(d, ReadCountMode.All);
				for (int v=0; v<r.getData().getVariationCount(d); v++) {
					if (r.getData().isMismatch(d, v)){//) && r.getData().getMismatchGenomic(d, v).charAt(0)=='T' && r.getData().getMismatchRead(d, v).charAt(0)=='C') {
						int mpos = r.getData().getMismatchPos(d, v);
						int pos = r.map(mpos);
						boolean overlap = r.getData().isPositionInOverlap(d, mpos);
						
						if (gene.getRegion().contains(pos) && (!overlap || counted.getQuick(mpos))) {
							counter.computeIfAbsent(pos, x->new double[2])[1]+=c;
						}
						counted.putQuick(mpos, true);
					}
				}
				
//				if (r.getReference().toString().equals("1-") && r.getRegion().contains(8963029))  
//					System.out.println(r+" "+c+" "+StringUtils.toString(counter.get(8963029)));
				
				cov.add(r.getRegion(), NumericArray.wrap(c));
				
			}
		}
		
//		if (gene.getReference().toString().equals("1-") && gene.getRegion().contains(8963029)) {
//			System.out.println(cov.getCoverages(gene.induce(8963029)));
//			System.out.println(cov.getProfile(0));
//		}
		
		senseCounter.addAndGet(co[0]);
		antisenseCounter.addAndGet(co[1]);
		
		
		StringBuilder sb = new StringBuilder();
		Iterator<Integer> it = counter.keySet().iterator();
		while (it.hasNext()) {
			int l = it.next();
			double[] count = counter.get(l);
//			if (l==25809577)
//				System.out.println(l);
			double c = cov.getCoverages(gene.induce(l)).getDouble(0);
			double pval = Beta.cumulative(conv, count[1]+conv, c-count[1]+1, true, false);
//			double l30 = Binomial.density(count[1], c, 0.3, true)/Math.log(10);
//			double l05 = Binomial.density(count[1], c, 0.05, true)/Math.log(10);
			
			// even if filtered much more stringently (likelihood of snp model must be 100x less than likelihood of conv model), the overall mismatch rates do not change (even for double)
			// regions with low coverage (like 1 or 2) might produce erroneous conversions (are in fact snps), i.e. bias the ntr towards 1; however low coverage genes are filtered anyways as not being expressed, and low coverage reagions in highly expressed genes contribute a negligible number of reads to the gene

			// this means that SNPs *MUST* be called with pooled data!
			
			
//			if (count[1]/c>cutoff && c>=coverageCutoff) { 
//			if (l05<2+l30) { 
			if (pval<pvalCutoff) {
				count[0] = c;
				
				sb.append(gene.getReference().toStrandIndependent()).append(":").append(l)
				.append("\t").append(String.format("%.1f", count[1]))
				.append("\t").append(String.format("%.1f", c))
				.append("\t").append(String.format("%.3g", pval))
				.append("\n");
				
			} else 
				it.remove();
			
			
			
		}
		if (sb.length()>0)
			sb.delete(sb.length()-1, sb.length());
		
		return sb.toString();
		
	}
	
	public static String collectNew(double pvalCutoff, GenomicRegionStorage<AlignedReadsData> reads, boolean[] no4sU, ImmutableReferenceGenomicRegion<String> gene, AtomicInteger senseCounter, AtomicInteger antisenseCounter) {
		gene=gene.toMutable().transformRegion(r->r.extendBack(1000).extendFront(1000)).toImmutable();
		
		
		CoverageAlgorithm cov = new CoverageAlgorithm(gene).setExample(NumericArray.wrap(0.0));
		
		int[] co = new int[2];
		ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> rit = 
				reads.ei(gene).sideEffect(r->co[0]++)
				.chain(reads.ei(gene.toMutable().toOppositeStrand()).sideEffect(r->co[1]++));
		BitVector counted = new BitVector(300);
		
		HashMap<Integer,double[]> counter = new HashMap<>();
		for (ImmutableReferenceGenomicRegion<AlignedReadsData> r : rit.loop()) {
			if (r.getRegion().getTotalLength()>counted.size()) 
				counted = new BitVector(r.getRegion().getTotalLength());
			
			for (int d=0; d<r.getData().getDistinctSequences(); d++) {
				counted.clear();
				double c = r.getData().getTotalCountForDistinct(d, ReadCountMode.All, (cond,count)->no4sU[cond]?count:0);
				
				for (int v=0; v<r.getData().getVariationCount(d); v++) {
					if (r.getData().isMismatch(d, v)){//) && r.getData().getMismatchGenomic(d, v).charAt(0)=='T' && r.getData().getMismatchRead(d, v).charAt(0)=='C') {
						int mpos = r.getData().getMismatchPos(d, v);
						int pos = r.map(mpos);
						boolean overlap = r.getData().isPositionInOverlap(d, mpos);
						
						if (gene.getRegion().contains(pos) && (!overlap || counted.getQuick(mpos))) {
							counter.computeIfAbsent(pos, x->new double[2])[1]+=c;
						}
						counted.putQuick(mpos, true);
					}
				}
				
//				if (r.getReference().toString().equals("1-") && r.getRegion().contains(8963029))  
//					System.out.println(r+" "+c+" "+StringUtils.toString(counter.get(8963029)));
				
				cov.add(r.getRegion(), NumericArray.wrap(c));
				
			}
		}
		
//		if (gene.getReference().toString().equals("1-") && gene.getRegion().contains(8963029)) {
//			System.out.println(cov.getCoverages(gene.induce(8963029)));
//			System.out.println(cov.getProfile(0));
//		}
		
		senseCounter.addAndGet(co[0]);
		antisenseCounter.addAndGet(co[1]);
		
		
		StringBuilder sb = new StringBuilder();
		Iterator<Integer> it = counter.keySet().iterator();
		while (it.hasNext()) {
			int l = it.next();
			double[] count = counter.get(l);

			int covcount = (int) cov.getCoverages(gene.induce(l)).getDouble(0);
			int convcount = (int) count[1];
			double pval = Binomial.cumulative(convcount-1, covcount, 0.001, false, false);

			if (pval<pvalCutoff) {
				count[0] = cov.getCoverages(gene.induce(l)).getDouble(0);
				
				sb.append(gene.getReference().toStrandIndependent()).append(":").append(l)
				.append("\t").append(String.format("%.1f", count[1]))
				.append("\t").append(String.format("%.1f", count[0]))
				.append("\t").append(String.format("%.3g", pval))
				.append("\n");
				
			} else 
				it.remove();
			
			
			
		}
		if (sb.length()>0)
			sb.delete(sb.length()-1, sb.length());
		
		return sb.toString();
		
	}
	
	public static void main(String[] args) {
		ImmutableReferenceGenomicRegion<String> gene = ImmutableReferenceGenomicRegion.parse("1+:1000-2000","Test");
		MemoryIntervalTreeStorage<DefaultAlignedReadsData> reads = new MemoryIntervalTreeStorage<>(DefaultAlignedReadsData.class);
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
		fac.start().newDistinctSequence().addMismatch(4, 'T', 'C', false).setCount(0, 1);
		reads.add(ImmutableReferenceGenomicRegion.parse("1+:1100-1110",fac.create()));
		
		CoverageAlgorithm cov = new CoverageAlgorithm(gene);
		
		HashMap<Integer,double[]> counter = new HashMap<>();
		for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> r : reads.ei(gene).loop()) {
			for (int d=0; d<r.getData().getDistinctSequences(); d++) {
				double c = r.getData().getTotalCountForDistinct(d, ReadCountMode.Weight);
				for (int v=0; v<r.getData().getVariationCount(d); v++) {
					if (r.getData().isMismatch(d, v) && r.getData().getMismatchGenomic(d, v).charAt(0)=='T' && r.getData().getMismatchRead(d, v).charAt(0)=='C') {
						int pos = r.map(r.getData().getMismatchPos(d, v));
						if (gene.getRegion().contains(pos))
							counter.computeIfAbsent(pos, x->new double[2])[1]+=c;
					}
				}
				cov.add(r.getRegion(), NumericArray.wrap(c));
			}
		}
		
		StringBuilder sb = new StringBuilder();
		Iterator<Integer> it = counter.keySet().iterator();
		while (it.hasNext()) {
			int l = it.next();
			double[] count = counter.get(l);
			
			double c = cov.getCoverages(gene.induce(l)).getDouble(0);
			
			if (count[1]/c>0.1) 
				count[0] = c;
			else 
				it.remove();
			
			sb.append(gene.getReference()).append(":").append(l)
				.append("\t").append(String.format("%.1f", count[0]))
				.append("\t").append(String.format("%.1f", count[1]))
				.append("\n");
			
		}
		System.out.println(sb);
	}
	
	
}
