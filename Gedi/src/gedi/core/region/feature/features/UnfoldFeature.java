package gedi.core.region.feature.features;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;

import java.util.Set;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class UnfoldFeature extends AbstractFeature<UnfoldGenomicRegionStatistics> {

	public UnfoldFeature() {
		minInputs = maxInputs = 1;
	}

	@Override
	protected void accept_internal(Set<UnfoldGenomicRegionStatistics> values) {
		Set<Object> in = getInput(0);
		values.add(()->in.iterator());
	}

	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		UnfoldFeature re = new UnfoldFeature();
		re.copyProperties(this);
		return re;
	}
	
}
