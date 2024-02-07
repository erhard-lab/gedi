package gedi.core.processing.old.transformer;

import gedi.core.processing.old.GenomicRegionProcessor;
import gedi.core.processing.old.ProcessorContext;
import gedi.core.region.MutableReferenceGenomicRegion;

public class RemoveIntronsProcessor implements GenomicRegionProcessor {


	
	@Override
	public void beginRegion(MutableReferenceGenomicRegion<?> region, ProcessorContext context)
			throws Exception {
		region.setRegion(region.getRegion().removeIntrons());
	}




}
