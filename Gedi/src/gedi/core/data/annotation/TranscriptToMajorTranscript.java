package gedi.core.data.annotation;

import java.util.Comparator;
import java.util.function.UnaryOperator;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;

public class TranscriptToMajorTranscript implements UnaryOperator<MemoryIntervalTreeStorage<Transcript>> {

	


	@Override
	public MemoryIntervalTreeStorage<Transcript> apply(
			MemoryIntervalTreeStorage<Transcript> storage) {
		
		Comparator<ImmutableReferenceGenomicRegion<Transcript>> cmp = (a,b)->a.getData().getGeneId().compareTo(b.getData().getGeneId());
		Comparator<ImmutableReferenceGenomicRegion<Transcript>> cmp2 = cmp.thenComparing((a,b)->Boolean.compare(b.getData().isCoding(), a.getData().isCoding()));
		cmp2 = cmp2.thenComparing((a,b)->Integer.compare(b.getRegion().getTotalLength(), a.getRegion().getTotalLength()));
		
		return storage.ei().
				sort(cmp2).
				multiplex(cmp, ImmutableReferenceGenomicRegion.class).
				map(a->a[0]).
				add(new MemoryIntervalTreeStorage<>(Transcript.class));
	}


}
