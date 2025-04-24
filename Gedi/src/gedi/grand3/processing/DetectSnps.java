package gedi.grand3.processing;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;
import java.util.logging.Logger;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.data.reads.subreads.MismatchReporterWithSequence;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.reads.MismatchPerPositionStatistics;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.TargetCollection;
import gedi.util.FunctorUtils.FilteredIterator;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.functions.CompoundParallelizedState;
import gedi.util.functions.ExtendedIterator;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.userInteraction.progress.Progress;
import jdistlib.Beta;

public class DetectSnps<A extends AlignedReadsData>  {

	private int nthreads;
	private ImmutableReferenceGenomicRegion<?> test;
	private double conv;
	private double pvalCutoff;
	
	private long[] senseAntisense = new long[3];
	
	private ReadSource<A> source;
	
	private boolean[] blackListConditions;

	
	public DetectSnps(ReadSource<A> source, double conv,
			double pvalCutoff, boolean[] blackListConditions) {
		this.source = source;
		this.conv = conv;
		this.pvalCutoff = pvalCutoff;
		this.blackListConditions = blackListConditions;
	}
	
	public DetectSnps<A>  setNthreads(int nthreads) {
		this.nthreads = nthreads;
		return this;
	}

	public void setTest(String loc) {
		this.test = ImmutableReferenceGenomicRegion.parse(loc);
	}

	public long getSense() {
		return senseAntisense[0];
	}
	public long getAntisense() {
		return senseAntisense[1];
	}
	
	/**
	 * The ReadSource should deliver reads from both strands
	 * @param progress
	 * @param targetRegions
	 * @param numTargetRegions
	 * @param source
	 * @param converter
	 * @param output
	 * @param posstat 
	 * @throws IOException
	 */
	public void process(Logger log, Supplier<Progress> progress, 
			TargetCollection targets,
			String output, 
			MismatchPerPositionStatistics posstat) throws IOException {
		
		AtomicLong senseCounter = new AtomicLong();
		AtomicLong antisenseCounter = new AtomicLong();
		CompoundParallelizedState state = new CompoundParallelizedState(posstat, NumericArray.createMemory(2, NumericArrayType.Integer));
		
		FilteredIterator<String> it = targets.iterateRegions()
			.iff(progress!=null,ei->ei.progress(progress.get(), targets.getNumRegions(), r->"Processing "+r.getData()))
			.parallelizedState(nthreads, 5, state, (ei,ss)->ei.map(l->{
				MismatchReporterWithSequence<MismatchPerPositionStatistics> pos = ((MismatchReporterWithSequence<MismatchPerPositionStatistics>) ss.get(0));
				NumericArray aa = ((NumericArray) ss.get(1));
				
				if (pos!=null)
					pos.cacheRegion(l, 100_000);
				
				
				int[] co = new int[2];
				ArrayList<ImmutableReferenceGenomicRegion<A>> list = source.getReads(l).sideEffect(r->{
					co[r.getReference().getStrand().equals(l.getReference().getStrand())?0:1]+=r.getData().getTotalCountOverallInt(ReadCountMode.Unique);
				}).list(); 
						
				// obtain all reads from both strands, in sense! does not matter, as mismatches per genomic position are counted, without paying attention to which mismatch it is!
				ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> mapped = 
						source.getConverter().convert(l.toLocationString()+":"+l.getData(),l.getReference(),list,pos,(u,t)->{aa.add(0, u); aa.add(1,t);});
				
				String re = findSnps(new ImmutableReferenceGenomicRegion<>(l.getReference().toStrandIndependent(),l.getRegion(),mapped));
				senseCounter.addAndGet(co[0]);
				antisenseCounter.addAndGet(co[1]);
				return re;
			}))
			.iff(progress!=null,ei->ei.progress(progress.get(), targets.getNumRegions(), r->"Finished"))
			.filter(s->s.length()>0);
		
		if (output!=null)
			it.print("Location\tCoverage\tMismatches\tP value",output);
		else
			it.drain();
		
		senseAntisense[0] = senseCounter.get();
		senseAntisense[1] = antisenseCounter.get();
		
		NumericArray aa = ((NumericArray) state.get(1));
		source.getConverter().logUsedTotal(log, aa.getInt(0), aa.getInt(1));
		
	}

	private String findSnps(ImmutableReferenceGenomicRegion<ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>>> cluster) {
		
		CoverageAlgorithm cov = new CoverageAlgorithm(cluster).setExample(NumericArray.wrap(0.0));
		
		HashMap<Integer,double[]> counter = new HashMap<>();
		HashSet<Integer> blacklisted = new HashSet<>();
		
		for (ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> r : cluster.getData().loop()) {
			for (int d=0; d<r.getData().getDistinctSequences(); d++) {
				double c = r.getData().getTotalCountForDistinct(d, ReadCountMode.All);
				double b = r.getData().getTotalCountForDistinct(d, ReadCountMode.All, (cond,count)->blackListConditions[cond]?count:0);

				for (int v=0; v<r.getData().getVariationCount(d); v++) {
					if (r.getData().isMismatch(d, v)){
						int mpos = r.getData().getMismatchPos(d, v);
						int pos = r.map(mpos);
						
						if (cluster.getRegion().contains(pos)) {
							if (b>0)
								blacklisted.add(pos);
							else
								counter.computeIfAbsent(pos, x->new double[2])[1]+=c;
						}
					}
				}
				cov.add(r.getRegion(), NumericArray.wrap(c));
			}
		}
		double[] BLACK = new double[0];
		for (Integer p : blacklisted)
			counter.put(p,BLACK);
		
		StringBuilder sb = new StringBuilder();
		Iterator<Integer> it = counter.keySet().iterator();
		while (it.hasNext()) {
			int l = it.next();
			
			double[] count = counter.get(l);
			
			if (count==BLACK) {
				sb.append(cluster.getReference()).append(":").append(l)
				.append("\tNA")
				.append("\tNA")
				.append("\t").append(String.format("%.3g", 0.0))
				.append("\n");
			}
			else {
				double c = cov.getCoverages(cluster.induce(l)).getDouble(0);
				double pval = Beta.cumulative(conv, count[1]+conv, c-count[1]+1, true, false);
				
	//			double l30 = Binomial.density(count[1], c, 0.3, true)/Math.log(10);
	//			double l05 = Binomial.density(count[1], c, 0.05, true)/Math.log(10);
				
				// even if filtered much more stringently (likelihood of snp model must be 100x less than likelihood of conv model), the overall mismatch rates do not change (even for double)
				// regions with low coverage (like 1 or 2) might produce erroneous conversions (are in fact snps), i.e. bias the ntr towards 1; however low coverage genes are filtered anyways as not being expressed, and low coverage reagions in highly expressed genes contribute a negligible number of reads to the gene
	
				// this means that SNPs *MUST* be called with pooled data!
				
				
	//			if (l05<2+l30) { 
				if (pval<pvalCutoff) {
					count[0] = c;
					
					sb.append(cluster.getReference()).append(":").append(l)
					.append("\t").append(String.format("%.1f", count[1]))
					.append("\t").append(String.format("%.1f", c))
					.append("\t").append(String.format("%.3g", pval))
					.append("\n");
					
				} //else 
					//it.remove();
			}
			
			
		}
		if (sb.length()>0)
			sb.delete(sb.length()-1, sb.length());
		
		return sb.toString();
		
	}

	
}
