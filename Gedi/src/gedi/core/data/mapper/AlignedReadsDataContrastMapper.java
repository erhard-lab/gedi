package gedi.core.data.mapper;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataMerger;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.ArrayUtils;
import gedi.util.FunctorUtils;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.mutable.MutableTriple;
import gedi.util.mutable.MutableTuple;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class AlignedReadsDataContrastMapper implements GenomicRegionDataMapper<IntervalTree<GenomicRegion,AlignedReadsData>, IntervalTree<GenomicRegion,AlignedReadsData>>{

	private ContrastMapping contrast;
	
	public AlignedReadsDataContrastMapper(ContrastMapping contrast) {
		this.contrast = contrast;
	}


	@Override
	public IntervalTree<GenomicRegion,AlignedReadsData> map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			IntervalTree<GenomicRegion, AlignedReadsData> data) {
		if (data.isEmpty()) return data;
		
		IntervalTree<GenomicRegion, AlignedReadsData> re = new IntervalTree<GenomicRegion,AlignedReadsData>(data.getReference());
		Iterator<Entry<GenomicRegion, AlignedReadsData>> it = data.entrySet().iterator();
		while (it.hasNext()) {
			Entry<GenomicRegion, AlignedReadsData> e = it.next();
			re.put(e.getKey(), new ConditionMappedAlignedReadsData(e.getValue(), contrast));
		}
		return re;
	}



}
