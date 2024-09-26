package gedi.core.data.reads;

import java.util.Iterator;
import java.util.function.Function;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.functions.EI;

public class UnfoldDistinct<R extends ReferenceGenomicRegion<? extends AlignedReadsData>> implements Function<R,Iterator<ImmutableReferenceGenomicRegion<SelectDistinctSequenceAlignedReadsData>>> {

	@Override
	public Iterator<ImmutableReferenceGenomicRegion<SelectDistinctSequenceAlignedReadsData>> apply(R t) {
		return EI.seq(0, t.getData().getDistinctSequences()).mapInt(d->new ImmutableReferenceGenomicRegion<SelectDistinctSequenceAlignedReadsData>(t.getReference(),t.getRegion(),new SelectDistinctSequenceAlignedReadsData(t.getData(),d)));
	}

}
