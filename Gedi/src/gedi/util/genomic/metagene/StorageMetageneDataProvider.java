package gedi.util.genomic.metagene;

import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.function.ToIntFunction;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

public class StorageMetageneDataProvider<D> implements MetageneDataProvider {

	private GenomicRegionStorage<D> storage;
	private Strandness strandness;
	private ToIntFunction<ReferenceGenomicRegion<D>> getPos;
	private ToIntFunction<ReferenceGenomicRegion<D>> getPos2;
	private ToDoubleFunction<ImmutableReferenceGenomicRegion<D>> getValue;
	private Predicate<ImmutableReferenceGenomicRegion<D>> filter;
	private Function<ImmutableReferenceGenomicRegion<?>,Predicate<ImmutableReferenceGenomicRegion<D>>> filterExt;
	
	/**
	 * 
	 * @param storage
	 * @param getPos returns coords in the genomic coordinate system (or -1 if invalid)
	 * @param getValue
	 */
	public StorageMetageneDataProvider(GenomicRegionStorage<D> storage, Strandness strandness, ToIntFunction<ReferenceGenomicRegion<D>> getPos,
			ToDoubleFunction<ImmutableReferenceGenomicRegion<D>> getValue) {
		this.storage = storage;
		this.strandness = strandness;
		this.getPos = getPos;
		this.getValue = getValue;
	}

	public StorageMetageneDataProvider(GenomicRegionStorage<D> storage, Strandness strandness, GenomicRegionPosition pos, int offset,
			ToDoubleFunction<ImmutableReferenceGenomicRegion<D>> getValue) {
		this.storage = storage;
		this.strandness = strandness;
		this.getPos = rgr->pos.position(rgr.getReference(), rgr.getRegion(), offset);
		this.getValue = getValue;
	}
	
	/**
	 * For coverage algorithm! must provide the first and last positions of the coverage
	 * @param storage
	 * @param strandness
	 * @param pos
	 * @param offset
	 * @param pos2
	 * @param offset2
	 * @param getValue
	 */
	public StorageMetageneDataProvider(GenomicRegionStorage<D> storage, Strandness strandness, GenomicRegionPosition pos, int offset,GenomicRegionPosition pos2, int offset2,
			ToDoubleFunction<ImmutableReferenceGenomicRegion<D>> getValue) {
		this.storage = storage;
		this.strandness = strandness;
		this.getPos = rgr->pos.position(rgr.getReference(), rgr.getRegion(), offset);
		this.getPos2 = rgr->pos2.position(rgr.getReference(), rgr.getRegion(), offset2);
		this.getValue = getValue;
	}
	
	
	public StorageMetageneDataProvider<D> setFilter(Predicate<ImmutableReferenceGenomicRegion<D>> filter) {
		this.filter = filter;
		return this;
	}
	
	public StorageMetageneDataProvider<D> setFilterExt(Function<ImmutableReferenceGenomicRegion<?>,Predicate<ImmutableReferenceGenomicRegion<D>>> filter) {
		this.filterExt = filter;
		return this;
	}
	
	
	@Override
	public AutoSparseDenseDoubleArrayCollector getData(ImmutableReferenceGenomicRegion<?> rgr) {
		AutoSparseDenseDoubleArrayCollector re = new AutoSparseDenseDoubleArrayCollector(rgr.getRegion().getTotalLength()/10,rgr.getRegion().getTotalLength());
		
		ExtendedIterator<ImmutableReferenceGenomicRegion<D>> it = EI.empty();
		if (strandness.equals(Strandness.Sense))
			it = storage.ei(rgr);
		else if (strandness.equals(Strandness.Antisense))
			it = storage.ei(rgr.toOppositeStrand());
		else {
			if (storage.getReferenceSequences().contains(rgr.getReference().toStrandIndependent()))
				it = storage.ei(rgr.toStrandIndependent());
			else
				it = storage.ei(rgr.toPlusStrand()).chain(storage.ei(rgr.toMinusStrand()));
		}
		
		if (filter!=null)
			it = it.filter(filter);
		if (filterExt!=null)
			it = it.filter(filterExt.apply(rgr));
		
		MutableReferenceGenomicRegion ee = rgr.toMutable();
		if (getPos2!=null) {
			for (ImmutableReferenceGenomicRegion<D> e : it.loop()) {
				ee.setRegion(e.getRegion());
				int pos = getPos.applyAsInt(ee);
				int pos2 = getPos2.applyAsInt(ee);
				if (pos>=0 && pos2>=0) {
					ArrayGenomicRegion fromTo = new ArrayGenomicRegion(pos,pos2+1).intersect(rgr.getRegion());
					if (!fromTo.isEmpty()) {
						fromTo = rgr.induce(fromTo);
						double value = getValue.applyAsDouble(e);
						for (int i=fromTo.getStart(); i<=fromTo.getStop(); i++)
							re.add(i, value);
						// difficult to speed up by the cumsum trick: if the stat is mean, the positive count might be added with length n, and the negative count subtracted with length n-1 (because the bins can have different lengths!)
					}
				}
			}
		}
		else {
			for (ImmutableReferenceGenomicRegion<D> e : it.loop()) {
				ee.setRegion(e.getRegion());
				int pos = getPos.applyAsInt(ee);
				if (pos>=0 && rgr.getRegion().contains(pos)) {
					pos = rgr.induce(pos);
					double value = getValue.applyAsDouble(e);
					re.add(pos, value);
				}
			}
		}
		return re;
	}

	
	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorageCoverage(GenomicRegionStorage<D> storage, Strandness strandness, 
			ReadCountMode mode, int... conditions) {
		return new StorageMetageneDataProvider<>(storage,strandness,GenomicRegionPosition.FivePrime,0,GenomicRegionPosition.ThreePrime,0,rgr->EI.wrap(conditions).mapToDouble(c->rgr.getData().getTotalCountForCondition(c, mode)).sum());
	}

	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorageCoverage(GenomicRegionStorage<D> storage, Strandness strandness, 
			ReadCountMode mode) {
		return new StorageMetageneDataProvider<>(storage,strandness,GenomicRegionPosition.FivePrime,0,GenomicRegionPosition.ThreePrime,0,rgr->rgr.getData().getTotalCountOverall(mode));
	}
	
	
	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorage(GenomicRegionStorage<D> storage, Strandness strandness, 
			GenomicRegionPosition pos, int offset,ReadCountMode mode, int... conditions) {
		return new StorageMetageneDataProvider<>(storage,strandness,pos,offset,rgr->EI.wrap(conditions).mapToDouble(c->rgr.getData().getTotalCountForCondition(c, mode)).sum());
	}

	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorage(GenomicRegionStorage<D> storage, Strandness strandness, 
			GenomicRegionPosition pos, int offset,ReadCountMode mode) {
		return new StorageMetageneDataProvider<>(storage,strandness,pos,offset,rgr->rgr.getData().getTotalCountOverall(mode));
	}
	
	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorage(GenomicRegionStorage<D> storage, Strandness strandness, 
			GenomicRegionPosition pos, ReadCountMode mode, int... conditions) {
		return new StorageMetageneDataProvider<>(storage,strandness,pos,0,rgr->EI.wrap(conditions).mapToDouble(c->rgr.getData().getTotalCountForCondition(c, mode)).sum());
	}

	public static <D extends AlignedReadsData> StorageMetageneDataProvider<D> fromReadStorage(GenomicRegionStorage<D> storage, Strandness strandness, 
			GenomicRegionPosition pos, ReadCountMode mode) {
		return new StorageMetageneDataProvider<>(storage,strandness,pos,0,rgr->rgr.getData().getTotalCountOverall(mode));
	}


	public static <D extends NumericArray> StorageMetageneDataProvider<D> fromCodonStorage(GenomicRegionStorage<D> storage, int frame, int condition) {
		frame = frame%3;
		if (frame<0) frame = 3+frame;
		int cframe = frame;
		StorageMetageneDataProvider re = new StorageMetageneDataProvider<>(storage,Strandness.Sense,GenomicRegionPosition.Start,0,cod->cod.getData().getDouble(condition));
		re.setFilterExt(orf->new Predicate<ImmutableReferenceGenomicRegion<NumericArray>>() {
			@Override
			public boolean test(ImmutableReferenceGenomicRegion<NumericArray> cod) {
				return ((ReferenceGenomicRegion) orf).induce(cod.getRegion()).getStart()%3==cframe;
			}
			
		});
		return re;
	}

	public static <D extends NumericArray> StorageMetageneDataProvider<D> fromCodonStorage(GenomicRegionStorage<D> storage, int condition) {
		return new StorageMetageneDataProvider<>(storage,Strandness.Sense,GenomicRegionPosition.FivePrime,0,cod->cod.getData().getDouble(condition));
	}

}
