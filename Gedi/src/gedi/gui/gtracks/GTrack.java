package gedi.gui.gtracks;

import gedi.core.data.mapper.DisablingGenomicRegionDataMapper;
import gedi.core.data.mapper.GenomicRegionDataMapper;
import gedi.gui.gtracks.rendering.GTracksRenderer;
import gedi.gui.gtracks.style.GTrackStyles;
import gedi.util.ReflectionUtils;
import gedi.util.dynamic.DynamicObject;

public interface GTrack<D,P> extends GenomicRegionDataMapper<D,GTracksRenderer>, DisablingGenomicRegionDataMapper<D, GTracksRenderer> {

	public enum HeightType {
		Fixed,// Will always render to a fixed height (e.g. sequence, packed tracks)
		Adaptive // Will always render into what is there
	}
	
	HeightType geHeightType();
	
	
	default void setStyles(DynamicObject styles) {
		try {
			ReflectionUtils.set(this, "style", new GTrackStyles(styles));
		} catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
		}
	}
	
}
