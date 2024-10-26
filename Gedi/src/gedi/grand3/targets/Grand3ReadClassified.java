package gedi.grand3.targets;

import java.util.Collection;
import java.util.List;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public interface Grand3ReadClassified {

	void classified(Collection<String> targets, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, CompatibilityCategory cat, ReadCountMode mode, boolean sense);
	
}
