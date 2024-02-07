package gedi.core.region.feature.features;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class ReadMismatchTypePositionFeature extends AbstractReadMismatchFeature {


	private String between1 = "->";
	private String between2 = " ";

	public void setBetween(String between1, String between2) {
		this.between1 = between1;
		this.between2 = between2;
	}


	@Override
	protected String getReturnValue(int pos, String g, String r) {
		if (g.equals(r))
			return g+between1+"N"+between2+pos;
		return g+between1+r+between2+pos;
	}

	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		ReadMismatchTypePositionFeature re = new ReadMismatchTypePositionFeature();
		re.copyProperties(this);
		re.between1 = between1;
		re.between2 = between2;
		return re;
	}

	

}
