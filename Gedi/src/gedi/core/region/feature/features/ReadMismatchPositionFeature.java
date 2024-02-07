package gedi.core.region.feature.features;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class ReadMismatchPositionFeature extends AbstractReadMismatchFeature {

	
	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		ReadMismatchPositionFeature re = new ReadMismatchPositionFeature();
		re.copyProperties(this);
		return re;
	}

	@Override
	protected String getReturnValue(int pos, String g, String r) {
		return pos+"";
	}

	

}
