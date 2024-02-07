package gedi.core.data.mapper;


import java.util.function.Consumer;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.dynamic.DynamicObject;

public interface GenomicRegionDataMapper<FROM,TO> {
	
	default void setJob(GenomicRegionDataMappingJob<FROM, TO> job){}
	default void setInput(int index, GenomicRegionDataMapper<?, FROM> input){}
	
	TO map(ReferenceSequence reference, GenomicRegion region, PixelLocationMapping pixelMapping, FROM data);
	default DynamicObject mapMeta(DynamicObject data) {
		return data;
	}
	
	default <T> void applyForAll(Class<T> cls, Consumer<T> consumer){}
	default boolean hasSideEffect() {
		return false;
	}
	
}
