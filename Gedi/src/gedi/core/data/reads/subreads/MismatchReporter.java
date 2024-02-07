package gedi.core.data.reads.subreads;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.IntToBooleanFunction;

public interface MismatchReporter {

	void startDistinct(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, int distinct, IntToBooleanFunction overlapperGetter);
	void reportMismatch(int variation, boolean overlap, boolean retained);
	default void reportConvertedReadsHistogram(int[][] readsHistoUsed, int[][] readsHistoDiscarded) {}

}
