package gedi.core.data.mapper;


import java.util.function.Consumer;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.dynamic.DynamicObject;

public interface GenomicRegionDataSource<TO> extends GenomicRegionDataMapper<Void,TO>  {
	
	
	TO get(ReferenceSequence reference, GenomicRegion region, PixelLocationMapping pixelMapping);
	default DynamicObject getMeta() {
		return DynamicObject.getEmpty();
	}

	default TO map(ReferenceSequence reference, GenomicRegion region,PixelLocationMapping pixelMapping,  Void from) {
		return get(reference,region,pixelMapping);
	}
	
	default DynamicObject mapMeta(DynamicObject meta) {
		return getMeta();
	}
	
	
	
}
