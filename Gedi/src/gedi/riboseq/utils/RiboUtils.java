package gedi.riboseq.utils;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.riboseq.inference.codon.Codon;
import gedi.util.FileUtils;
import gedi.util.ParseUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.tree.redblacktree.IntervalTreeSet;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.parsing.ReferenceGenomicRegionParser;
import gedi.util.userInteraction.progress.Progress;

import java.io.IOException;
import java.util.HashSet;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;

import cern.colt.bitvector.BitVector;

public class RiboUtils {

	public static boolean isLeadingMismatchInsideGenomicRegion(AlignedReadsData ard, int distinct) {
		for (int i=0; i<ard.getVariationCount(distinct); i++) {
			if (ard.isMismatch(distinct, i) && ard.getMismatchPos(distinct, i)==0)
				return true;
			if (ard.isSoftclip(distinct, i) && ard.isSoftclip5p(distinct, i) && ard.getSoftclip(distinct, i).length()==1)
				return false;
		}
		throw new RuntimeException("No leading mismatch!");
	}
	
	public static boolean hasLeadingMismatch(AlignedReadsData ard, int distinct) {
		for (int i=0; i<ard.getVariationCount(distinct); i++) {
			if (ard.isMismatch(distinct, i) && ard.getMismatchPos(distinct, i)==0)
				return true;
			if (ard.isSoftclip(distinct, i) && ard.isSoftclip5p(distinct, i) && ard.getSoftclip(distinct, i).length()==1)
				return true;
		}
		return false;
	}
	
	public static char getLeadingMismatch(AlignedReadsData ard, int distinct) {
		for (int i=0; i<ard.getVariationCount(distinct); i++) {
			if (ard.isMismatch(distinct, i) && ard.getMismatchPos(distinct, i)==0)
				return ard.getMismatchRead(distinct, i).charAt(0);
			if (ard.isSoftclip(distinct, i) && ard.isSoftclip5p(distinct, i) && ard.getSoftclip(distinct, i).length()==1)
				return 'N';
		}
		return '\0';
	}

	public static Predicate<ReferenceGenomicRegion<AlignedReadsData>> parseReadFilter(String spec) {
		if (spec==null||spec.length()==0) return null;
		
		if (spec.startsWith("[") && spec.endsWith("]")) {
			HashSet<String> refs = EI.wrap(StringUtils.split(spec.substring(1, spec.length()-1), ',')).map(s->StringUtils.trim(s)).toCollection(new HashSet<String>());
			return rgr->refs.contains(rgr.getReference().getName());
		}
		ReferenceGenomicRegionParser<Void> p = new ReferenceGenomicRegionParser<Void>();
		if (p.canParse(spec)) {
			MutableReferenceGenomicRegion ref = p.apply(spec);
			if (ref.getRegion().getTotalLength()>1)
				return read->ref.intersects(read);
		}
		BitVector bv = ParseUtils.parseRangeBv(spec, 500);
		return rgr->bv.getQuick(rgr.getRegion().getTotalLength());
	}

	public static boolean isInframeOverlap(GenomicRegion a, GenomicRegion b) {
		ArrayGenomicRegion o = a.intersect(b);
		if (o.isEmpty()) return false;
		return a.induce(o).getStart()%3==b.induce(o).getStart()%3 || (a.getTotalLength()-a.induce(o).getEnd())%3==(b.getTotalLength()-b.induce(o).getEnd())%3; 
	}
	
	
	public static <A> void processCodonsSink(String file, String title, Progress progress, Supplier<String> extraProgress, int nthreads, int chunk, Class<A> aclass, 
			Function<ExtendedIterator<ImmutableReferenceGenomicRegion<IntervalTreeSet<Codon>>>,ExtendedIterator<ImmutableReferenceGenomicRegion<A>>> parallel, 
			Consumer<ImmutableReferenceGenomicRegion<A>> sink) throws IOException {
	
		processCodons(file, title, progress, extraProgress, nthreads, chunk, aclass, parallel, oit->oit.forEachRemaining(sink));
		
	}
	
	public static <A> void processCodons(String file, String title, Progress progress, Supplier<String> extraProgress, int nthreads, int chunk, Class<A> aclass, 
			Function<ExtendedIterator<ImmutableReferenceGenomicRegion<IntervalTreeSet<Codon>>>,ExtendedIterator<ImmutableReferenceGenomicRegion<A>>> parallel, 
			Consumer<ExtendedIterator<ImmutableReferenceGenomicRegion<A>>> sink) throws IOException {
		
		PageFile prof = new PageFile(file);
		Progress uprog = progress;
		progress.init();
		ExtendedIterator<ImmutableReferenceGenomicRegion<A>> oit = prof.ei().map(pr->{
					try {
						ReferenceSequence ref = FileUtils.readReferenceSequence(pr);
						GenomicRegion reg = FileUtils.readGenomicRegion(pr);
						int n = pr.getCInt();
						IntervalTreeSet<Codon> codons = new IntervalTreeSet<Codon>(null);
						for (int i=0; i<n;i++) {
							codons.add(new Codon(FileUtils.readGenomicRegion(pr),FileUtils.readDoubleArray(pr)));
						}
						
						return new ImmutableReferenceGenomicRegion<IntervalTreeSet<Codon>>(ref,reg,codons);
					} catch (IOException e) {
						throw new RuntimeException("Could not read temporary profiles!",e);
					}
					})
		.parallelized(nthreads, chunk, ei->parallel.apply(ei)
				.sideEffect(r->{
					synchronized (uprog) {
						
						
						if (extraProgress!=null)
							uprog.setDescription(()->title+" "+r.toLocationStringRemovedIntrons()+" mem="+StringUtils.getShortHumanReadableMemory(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())
									+" "+extraProgress.get()).incrementProgress();
						else
							uprog.setDescription(()->title+" "+r.toLocationStringRemovedIntrons()).incrementProgress();
					}
				})
				);
		
		sink.accept(oit);
		
		progress.finish();
		
	}

	
}
