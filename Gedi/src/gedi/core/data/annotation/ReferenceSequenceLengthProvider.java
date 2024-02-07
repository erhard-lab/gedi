package gedi.core.data.annotation;

import gedi.core.region.ArrayGenomicRegion;

public interface ReferenceSequenceLengthProvider {

	/**
	 * Negative numbers: length is lower bound of true length
	 * @param reference
	 * @return
	 */
	int getLength(String name);

	default ArrayGenomicRegion getRegionOfReference(String name) {
		return new ArrayGenomicRegion(0,getLength(name));
	}
	
}
