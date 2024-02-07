package gedi.core.data.mapper;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;

public interface DisablingGenomicRegionDataMapper<FROM,TO> extends GenomicRegionDataMapper<FROM, TO> {
	
	boolean isDisabled(ReferenceSequence reference, GenomicRegion region, PixelLocationMapping pixelMapping);

	
}
