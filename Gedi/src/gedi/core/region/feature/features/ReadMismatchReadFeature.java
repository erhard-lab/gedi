package gedi.core.region.feature.features;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class ReadMismatchReadFeature extends AbstractReadMismatchFeature {



	@Override
	protected String getReturnValue(int pos, String g, String r) {
		if (g.equals(r))
			return "N";
		return r;
	}

	
	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		ReadMismatchReadFeature re = new ReadMismatchReadFeature();
		re.copyProperties(this);
		return re;
	}

	

}
