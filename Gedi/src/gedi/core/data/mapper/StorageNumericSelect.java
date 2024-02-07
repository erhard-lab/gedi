package gedi.core.data.mapper;

import gedi.util.datastructure.tree.redblacktree.IntervalTree;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class StorageNumericSelect extends StorageNumericCompute{


	public StorageNumericSelect(String range) {
		super(new SelectOp(range));
	}

	
}
