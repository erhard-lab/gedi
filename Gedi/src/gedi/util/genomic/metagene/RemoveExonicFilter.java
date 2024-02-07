package gedi.util.genomic.metagene;

import java.util.function.Predicate;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;

public class RemoveExonicFilter implements Predicate<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>>{

	private MemoryIntervalTreeStorage<Transcript> transcripts;
	
	
	public RemoveExonicFilter(MemoryIntervalTreeStorage<Transcript> transcripts) {
		this.transcripts = transcripts;
	}
	public RemoveExonicFilter(Genomic g) {
		this.transcripts = g.getTranscripts();
	}
	

	@Override
	public boolean test(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> r) {
		if (r.getData().getNumParts(r, 0)>1) return false;
		int maxOverlap = transcripts.ei(r).mapToInt(t->t.getRegion().intersectLength(r.getRegion())).max();
		return maxOverlap<r.getRegion().getTotalLength()-3;
	}

}
