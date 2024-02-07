package gedi.core.data.mapper;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.ParseUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class AlignedReadsDataSelectMapper implements GenomicRegionDataMapper<IntervalTree<GenomicRegion, AlignedReadsData>, IntervalTree<GenomicRegion, AlignedReadsData>> {

    private String range;

    public AlignedReadsDataSelectMapper(String range) {
        this.range = range;
    }

    @Override
    public IntervalTree<GenomicRegion, AlignedReadsData> map(ReferenceSequence reference, GenomicRegion region, PixelLocationMapping pixelMapping, IntervalTree<GenomicRegion, AlignedReadsData> data) {
        int numCond = data.ei().next().getData().getNumConditions();

        int[] ranges = ParseUtils.parseRangePositions(range, numCond, new IntArrayList(), null, ':', false).toIntArray();
        ContrastMapping mapping = new ContrastMapping();
        for (int i=0; i<ranges.length; i++)
            mapping.addMapping(ranges[i], i);


        IntervalTree<GenomicRegion,AlignedReadsData> re = new IntervalTree<>(data.getReference());
        data.ei().forEachRemaining(rgr -> re.put(rgr.getRegion(), new ConditionMappedAlignedReadsData(rgr.getData(), mapping)));

        return re;
    }
}
