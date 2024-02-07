package gedi.core.data.mapper;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataMerger;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.ArrayUtils;
import gedi.util.FunctorUtils;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.functions.EI;
import gedi.util.mutable.MutableTuple;

@GenomicRegionDataMapping(fromType=MutableTuple.class,toType=IntervalTree.class)
public class AlignedReadsDataMergeMapper implements GenomicRegionDataMapper<MutableTuple, IntervalTree<GenomicRegion,AlignedReadsData>>{

	private ArrayList<int[]> mapping = new ArrayList<>();
	private int[] numConditions;
	private ContrastMapping contrast;
	
	
	
	public AlignedReadsDataMergeMapper(int[] numConditions) {
		this.numConditions = numConditions;
	}
	
	
	public void map(int file, int condition, int name) {
		mapping.add(new int[] {file,condition,name});
	}

	
	@Override
	public IntervalTree<GenomicRegion,AlignedReadsData> map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			MutableTuple data) {
		
		if (contrast==null){
			synchronized (this) {
				if (contrast==null) {
					contrast = new ContrastMapping();
					int[] cumNumConditions = ArrayUtils.cumSum(numConditions, 1);
					for (int[] tr : mapping) {
						contrast.addMapping((tr[0]>0?cumNumConditions[tr[0]-1]:0)+tr[1], tr[2]);
					}
					contrast.build();
				}
			}
		}
		
		Iterator<ImmutableReferenceGenomicRegion<AlignedReadsData>>[] its = new Iterator[data.size()];
		for (int i=0; i<its.length; i++) {
			IntervalTree<GenomicRegion, AlignedReadsData> itree = data.<IntervalTree<GenomicRegion,AlignedReadsData>>get(i);
			if (itree==null || itree.isEmpty()) {
				its[i] = EI.empty();
			} else {
				its[i] = itree.ei();
			}
		}
		
		
		if (its.length>1) {
			AlignedReadsDataMerger merger = new AlignedReadsDataMerger(numConditions);
			IntervalTree<GenomicRegion,AlignedReadsData> re = new IntervalTree<>(data.<IntervalTree<GenomicRegion,AlignedReadsData>>get(0).getReference());
			
			Comparator<ImmutableReferenceGenomicRegion<AlignedReadsData>> comp = FunctorUtils.<ImmutableReferenceGenomicRegion<AlignedReadsData>>naturalComparator();
			Class<ImmutableReferenceGenomicRegion<AlignedReadsData>> cls = (Class)ImmutableReferenceGenomicRegion.class;
			for (ImmutableReferenceGenomicRegion<AlignedReadsData>[] a : FunctorUtils.parallellIterator(its, comp, cls).loop()) {
				ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> r = merger.merge(a);
				ConditionMappedAlignedReadsData mdata = new ConditionMappedAlignedReadsData(r.getData(), contrast);
				if (mdata.getTotalCountOverall(ReadCountMode.All)>0)
					re.put(r.getRegion(), mdata);
			}
			return re;
			
		} else {
			IntervalTree<GenomicRegion,AlignedReadsData> re = new IntervalTree<>(data.<IntervalTree<GenomicRegion,AlignedReadsData>>get(0).getReference());
			while (its[0].hasNext()) {
				ImmutableReferenceGenomicRegion<AlignedReadsData> r = its[0].next();
				ConditionMappedAlignedReadsData mdata = new ConditionMappedAlignedReadsData(r.getData(), contrast);
				if (mdata.getTotalCountOverall(ReadCountMode.All)>0)
					re.put(r.getRegion(), mdata);
			}
			
			return re;
		}
			
		
		
	}



}
