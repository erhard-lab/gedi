package gedi.core.data.mapper;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.SelectDistinctSequenceAlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class AlignedReadsFilter extends StorageFilter<AlignedReadsData>{

	

	public void unique() {
		addFilter(r->r.getData().isAnyUniqueMapping());
		addOperator(r->{
			if (!r.getData().isAnyAmbigousMapping()) return r;
			IntArrayList ds = new IntArrayList(r.getData().getDistinctSequences());
			for (int d=0; d<r.getData().getDistinctSequences(); d++)
				if (r.getData().isUniqueMapping(d))
					ds.add(d);
			ReferenceGenomicRegion<AlignedReadsData> re = new ImmutableReferenceGenomicRegion<>(r.getReference(),r.getRegion(),new SelectDistinctSequenceAlignedReadsData(r.getData(), ds.toIntArray()));
			return re;
		});
	}
	

}
