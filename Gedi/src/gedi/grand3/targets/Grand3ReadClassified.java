package gedi.grand3.targets;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public interface Grand3ReadClassified {

	void classified(ImmutableReferenceGenomicRegion<String> target, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, CompatibilityCategory cat, ReadCountMode mode, boolean sense);
	
}
