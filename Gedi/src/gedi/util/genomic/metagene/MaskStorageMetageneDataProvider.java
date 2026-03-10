package gedi.util.genomic.metagene;

import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.EI;

public class MaskStorageMetageneDataProvider<D> implements MetageneDataProvider {

	private MemoryIntervalTreeStorage<?> storage;
	

	public MaskStorageMetageneDataProvider(ImmutableReferenceGenomicRegion<?>[] mask) {
		this(mask,0,0);
	}
	public MaskStorageMetageneDataProvider(ImmutableReferenceGenomicRegion<?>[] mask, int offsetUpstream, int offsetDownstream) {
		storage = new MemoryIntervalTreeStorage(Void.class);
		for (ImmutableReferenceGenomicRegion r : mask) {
			if (offsetUpstream!=0 || offsetDownstream!=0) {
				if (r.getRegion().getNumParts()>1) throw new RuntimeException("Not defined, extending and introns in mask!");
				GenomicRegion reg = r.getReference().isPlus()?
						new ArrayGenomicRegion(r.getRegion().getStart()+offsetUpstream,r.getRegion().getEnd()+offsetDownstream):
							new ArrayGenomicRegion(r.getRegion().getStart()-offsetDownstream,r.getRegion().getEnd()-offsetUpstream); 
				r = new ImmutableReferenceGenomicRegion<D>(r.getReference(), reg);
			}
			storage.add(r);
		}
	}
	
	@Override
	public AutoSparseDenseDoubleArrayCollector getData(ImmutableReferenceGenomicRegion<?> rgr) {
		AutoSparseDenseDoubleArrayCollector re = new AutoSparseDenseDoubleArrayCollector(rgr.getRegion().getTotalLength()/10,rgr.getRegion().getTotalLength());
		
		for (ImmutableReferenceGenomicRegion<?> tr : storage.ei(rgr).loop()) {
			for (int p=0; p<tr.getRegion().getNumParts(); p++) {
				for (int i=tr.getRegion().getStart(p); i<tr.getRegion().getEnd(p); i++)
					if (rgr.getRegion().contains(i))
						re.add(rgr.induce(i),1);
			}
		}

		return re;
	}

	
	
}
