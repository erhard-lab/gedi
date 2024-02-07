package gedi.core.region.feature.features;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.special.UnfoldGenomicRegionStatistics;
import gedi.util.functions.EI;

import java.util.ArrayList;
import java.util.Set;


@GenomicRegionFeatureDescription(toType=UnfoldGenomicRegionStatistics.class)
public class ReadMismatchGenomicFeature extends AbstractReadMismatchFeature {


	@Override
	protected String getReturnValue(int pos, String g, String r) {
		return g;
	}


	@Override
	public GenomicRegionFeature<UnfoldGenomicRegionStatistics> copy() {
		ReadMismatchGenomicFeature re = new ReadMismatchGenomicFeature();
		re.copyProperties(this);
		return re;
	}

	

}
