package gedi.core.region.feature.features;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;

import java.util.Set;


@GenomicRegionFeatureDescription(fromType=ReferenceGenomicRegion.class,toType=ReferenceGenomicRegion.class)
public class ContainedFeature extends AbstractFeature<ReferenceGenomicRegion<?>> {

	
	
	
	
	
	public ContainedFeature() {
		minInputs = maxInputs = 1;
	}
	
	@Override
	public GenomicRegionFeature<ReferenceGenomicRegion<?>> copy() {
		ContainedFeature re = new ContainedFeature();
		re.copyProperties(this);
		return re;
	}


	@Override
	protected void accept_internal(Set<ReferenceGenomicRegion<?>> values) {
	
		Set<ReferenceGenomicRegion<?>> inputs = getInput(0);
		
		for (ReferenceGenomicRegion<?> rgr : inputs) {
			if (referenceRegion.getData() instanceof AlignedReadsData) {
				AlignedReadsData d = (AlignedReadsData)referenceRegion.getData();
				if (d.isConsistentlyContained(referenceRegion,rgr,0))
					values.add(rgr);
			}
			else if (rgr.getRegion().containsUnspliced(referenceRegion.getRegion()))
				values.add(rgr);
			
		}
		
	}
	

}
