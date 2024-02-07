package gedi.core.data.annotation;

import gedi.core.region.ReferenceGenomicRegion;

/**
 * If a data object of a {@link ReferenceGenomicRegion} has position information, this must be mapped along with the region!
 * @author erhard
 *
 */
public interface GenomicRegionMappable<T extends GenomicRegionMappable<T>> {

	T map(ReferenceGenomicRegion<?> mapper);
	T induce(ReferenceGenomicRegion<?> mapper);
	
}
