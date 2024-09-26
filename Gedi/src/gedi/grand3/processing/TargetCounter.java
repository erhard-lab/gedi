package gedi.grand3.processing;

import java.util.List;

import gedi.core.region.ImmutableReferenceGenomicRegion;

public interface TargetCounter<T extends TargetCounter<T,R>,R> extends SubreadCounter<T> {

	void startTarget(String currentTarget);
	List<R> getResultForCurrentTarget();
	
	boolean mergeWithSameName();
	
}
