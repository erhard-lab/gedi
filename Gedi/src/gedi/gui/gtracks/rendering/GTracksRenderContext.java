package gedi.gui.gtracks.rendering;

import java.util.LinkedHashMap;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.gui.PixelBasepairMapper;

public class GTracksRenderContext {

	
	private LinkedHashMap<ReferenceSequence, GenomicRegion> regions = new LinkedHashMap<>();
	private PixelBasepairMapper xmapper;

	public GTracksRenderContext(ReferenceSequence[] references, GenomicRegion[] regions, PixelBasepairMapper xmapper) {
		for (int i=0; i<references.length; i++)
			this.regions.put(references[i],regions[i]);
		this.xmapper = xmapper;
	}
	
	public GTracksRenderContext(ReferenceSequence reference, GenomicRegion region, PixelBasepairMapper xmapper) {
		this.regions.put(reference,region);
		this.xmapper = xmapper;
	}
	
	
	public GenomicRegion getRegionToRender(ReferenceSequence reference) {
		return regions.get(reference);
	}

	public PixelBasepairMapper getLocationMapper() {
		return xmapper;
	}

}
