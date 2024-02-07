package gedi.grand3.processing;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.targets.TargetCollection;
import gedi.util.functions.ParallelizedState;

public interface SubreadCounter<T extends SubreadCounter<T>> extends ParallelizedState<T> {

	/**
	 * This is called per distinct sequence (only if {@link TargetCollection#classify(ImmutableReferenceGenomicRegion, ImmutableReferenceGenomicRegion, boolean, gedi.grand3.targets.Grand3ReadClassified)} is ok with it!
	 * @param buffer
	 */
	void count(SubreadProcessorMismatchBuffer buffer);

}
