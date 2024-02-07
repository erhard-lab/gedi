package gedi.core.region.feature.features;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class ReadMismatchTypeFeature extends AbstractReadMismatchFeature {


	private String between = "->";

	public void setBetween(String between) {
		this.between = between;
	}


	@Override
	protected String getReturnValue(int pos, String g, String r) {
		if (g.equals(r))
			return g+between+"N";
		return g+between+r;
	}

	
	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		ReadMismatchTypeFeature re = new ReadMismatchTypeFeature();
		re.copyProperties(this);
		re.between = between;
		return re;
	}

	

}
