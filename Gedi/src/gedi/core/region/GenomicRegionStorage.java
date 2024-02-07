package gedi.core.region;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.function.BiFunction;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;

import gedi.core.data.HasConditions;
import gedi.core.data.annotation.ReferenceSequencesProvider;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoaderExtensionPoint;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.MappedSpliterator;
import gedi.util.functions.SpliteratorArraySpliterator;
import gedi.util.userInteraction.progress.Progress;


public interface GenomicRegionStorage<D> extends ReferenceSequencesProvider {

	public static <D> Spliterator<ImmutableReferenceGenomicRegion<D>> immut(Spliterator<MutableReferenceGenomicRegion<D>> split) {
		return new MappedSpliterator<MutableReferenceGenomicRegion<D>,ImmutableReferenceGenomicRegion<D>>(split,m->m.toImmutable());
	}
	
	public static <D> Spliterator<GenomicRegion> toRegion(Spliterator<MutableReferenceGenomicRegion<D>> split) {
		return new MappedSpliterator<MutableReferenceGenomicRegion<D>,GenomicRegion>(split,m->m.getRegion());
	}
	
	
	
	/**
	 * Available slots:
	 * 
	 * for AlignedReadsData storages:
	 * .conditions: Array of length {@link AlignedReadsData#getNumConditions()}, slots: name,color
	 * 
	 * 
	 * @return
	 */
	DynamicObject getMetaData();
	void setMetaData(DynamicObject meta);
	
	/**
	 * Gets all {@link ReferenceSequence} contained in this storage. Note that there is not necessarily any region associated with 
	 * each returned reference!
	 * @return
	 */
	Set<ReferenceSequence> getReferenceSequences();
	
	default ExtendedIterator<ReferenceSequence> iterateReferenceSequences() {
		return EI.wrap(getReferenceSequences());
	}
	
	Class<D> getType();
	
	/**
	 * Iterates over all {@link ImmutableReferenceGenomicRegion}s in this storage. Does not return null, but an empty spliterator, when applicable.
	 * @return
	 */
	default Spliterator<MutableReferenceGenomicRegion<D>> iterateMutableReferenceGenomicRegions() {
		ReferenceSequence[] refs = getReferenceSequences().toArray(new ReferenceSequence[0]);
		Supplier<Spliterator<MutableReferenceGenomicRegion<D>>>[] re = new Supplier[refs.length];
		Arrays.sort(refs);
		for (int i=0; i<re.length; i++) {
			int index = i;
			re[i] = ()->iterateMutableReferenceGenomicRegions(refs[index]);
		}
		return new SpliteratorArraySpliterator<MutableReferenceGenomicRegion<D>>(re);
	}
	
	default MemoryIntervalTreeStorage<D> toMemory() {
		MemoryIntervalTreeStorage<D> re = new MemoryIntervalTreeStorage<D>(getType());
		re.setMetaData(getMetaData());
		re.fill(this);
		return re;
	}
	
	
	/**
	 * Iterates over all {@link ImmutableReferenceGenomicRegion}s belonging to the given reference in this storage. Does not return null, but an empty spliterator, when applicable.
	 * @return
	 */
	Spliterator<MutableReferenceGenomicRegion<D>> iterateMutableReferenceGenomicRegions(ReferenceSequence reference);

	/**
	 * Iterates over all {@link ImmutableReferenceGenomicRegion}s belonging to the given reference and intersecting the given region in this storage. Does not return null, but an empty spliterator, when applicable.
	 * Two regions intersect, if any of the parts intersect.
	 * 
	 * @param reference
	 * @param region
	 * @return
	 */
	Spliterator<MutableReferenceGenomicRegion<D>> iterateIntersectingMutableReferenceGenomicRegions(ReferenceSequence reference, GenomicRegion region);
	
	/**
	 * Adds the given region. Note that the storage may not support adding elements or may only support adding elements in a specific order!
	 * If the region already is in this storage, the data gets overwritten 
	 * 
	 * @see #add(ReferenceSequence, GenomicRegion, Object, BinaryOperator)
	 * @param reference
	 * @param region
	 * @param data
	 * @return true if the region has not been in the storage prior to add, and false, if
	 * 			it already was in it (irrespective of which data was associated with it).
	 */
	boolean add(ReferenceSequence reference, GenomicRegion region, D data);
	
	/**
	 * Can be a file, a folder, or some connection string.
	 * Only null, if the data is in main memory!
	 * @return
	 */
	String getPath();
	
	default boolean add(ReferenceGenomicRegion<D> rgr) {
		return add(rgr.getReference(),rgr.getRegion(),rgr.getData());
	}
	
	/**
	 * Removes the given region. Note that the storage may not support removing elements!
	 * @param reference
	 * @param region
	 * @param data
	 * @return true if the region has been in the storage prior to remove.
	 */
	boolean remove(ReferenceSequence reference, GenomicRegion region);
	
	/**
	 * Gets whether the region is in the storage. Note that this is not the same as {@link #getData(ReferenceSequence, GenomicRegion)}!=null,
	 * since null values are valid data entries for existing regions.
	 * @param reference
	 * @param region
	 * @return
	 */
	boolean contains(ReferenceSequence reference, GenomicRegion region);
	
	/**
	 * Gets the data associated with the given region. Note that the return value null can mean two things: Either the region is not in this storage,
	 * or it is in it with null associated with it.
	 * @param reference
	 * @param region
	 * @return
	 */
	D getData(ReferenceSequence reference, GenomicRegion region);
	
	long size(ReferenceSequence reference);
	void clear();
	
	
	default long size() {
		long re = 0;
		for (ReferenceSequence ref : getReferenceSequences())
			re+=size(ref);
		return re;
	}
	
	default void fill(Iterator<? extends ReferenceGenomicRegion<D>> it) {
		fill(it,null);
	}

	default void fill(GenomicRegionStorage<D> storage) {
		fill(storage,(Progress)null);
	}
	
	default <O> void fill(GenomicRegionStorage<O> storage, Function<MutableReferenceGenomicRegion<O>,MutableReferenceGenomicRegion<D>> mapper) {
		fill(storage,mapper,null);
	}
	
	default void fill(Iterator<? extends ReferenceGenomicRegion<D>> it, Progress p) {
		while (it.hasNext()) {
			ReferenceGenomicRegion<D> rgr = it.next();
			add(rgr.getReference(),rgr.getRegion(),rgr.getData());
		}
	}

	default void fill(GenomicRegionStorage<D> storage, Progress p) {
		fill(storage,t->t,p);
	}
	
	default <O> void fill(GenomicRegionStorage<O> storage, Function<MutableReferenceGenomicRegion<O>,MutableReferenceGenomicRegion<D>> mapper, Progress p) {
		fill(EI.wrap(Spliterators.iterator(storage.iterateMutableReferenceGenomicRegions())).map(mapper),null);
	}

	/**
	 * Adds the given region. Note that the storage may not support adding elements or may only support adding elements in a specific order!
	 * If the region already is in this storage, the data is combined using the given combiner  
	 * 
	 * @param reference
	 * @param region
	 * @param data
	 * @param combiner
	 * @return true if the region has not been in the storage prior to add, and false, if
	 * 			it already was in it (irrespective of which data was associated with it).
	 */
	default boolean add(ReferenceSequence reference, GenomicRegion region, D data, BinaryOperator<D> combiner) {
		if (contains(reference, region)) 
			return add(reference,region,combiner.apply(getData(reference, region), data));
		return add(reference,region,data);
	}
	
	
	default Spliterator<MutableReferenceGenomicRegion<D>> iterateIntersectingMutableReferenceGenomicRegions(ReferenceSequence reference, int start, int end) {
		return iterateIntersectingMutableReferenceGenomicRegions(reference, new ArrayGenomicRegion(start,end));
	}
	
	default Spliterator<ImmutableReferenceGenomicRegion<D>> iterateReferenceGenomicRegions() {
		return immut(iterateMutableReferenceGenomicRegions());
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei() {
		return EI.wrap(iterateReferenceGenomicRegions());
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(ReferenceSequence reference) {
		return EI.wrap(iterateReferenceGenomicRegions(reference));
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(ReferenceSequence... reference) {
		if (reference.length==0) return EI.empty();
		ExtendedIterator<ImmutableReferenceGenomicRegion<D>> re = EI.wrap(iterateReferenceGenomicRegions(reference[0]));
		for (int i=1; i<reference.length; i++) 
			re = re.chain(EI.wrap(iterateReferenceGenomicRegions(reference[i])));
		return re;
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(ReferenceSequence reference, GenomicRegion region) {
		return EI.wrap(immut(iterateIntersectingMutableReferenceGenomicRegions(reference, region)));
	}
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(ReferenceGenomicRegion<?> rgr) {
		if (rgr==null)
			return ei();
		return EI.wrap(immut(iterateIntersectingMutableReferenceGenomicRegions(rgr.getReference(), rgr.getRegion())));
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(String location) {
		if (location==null)
			return ei();
		location = StringUtils.trim(location);
		
		if (location.contains(";"))
			return EI.split(location, ';').unfold(loc->ei(loc));
		if (location.contains(","))
			return EI.split(location, ',').unfold(loc->ei(loc));
		
		if (!location.contains(":")) {
			Chromosome ref = Chromosome.obtain(location);
			if (getReferenceSequences().contains(ref) || !ref.getStrand().equals(Strand.Independent))
				return ei(ref);
			return ei(ref.toPlusStrand()).chain(ei(ref.toMinusStrand()));
		}
		return ei(ImmutableReferenceGenomicRegion.parse(location));
	}
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<D>> ei(String...location) {
		if (location.length==0) return EI.empty();
		ExtendedIterator<ImmutableReferenceGenomicRegion<D>> re = ei(location[0]);
		for (int i=1; i<location.length; i++) 
			re = re.chain(ei(location[i]));
		return re;
	}
	
	default  <C extends Collection<ImmutableReferenceGenomicRegion<D>>> C getReferenceGenomicRegions(C re) {
		iterateMutableReferenceGenomicRegions().forEachRemaining(mrgr->re.add(mrgr.toImmutable()));
		return re;
	}
	
	default ArrayList<ImmutableReferenceGenomicRegion<D>> getReferenceGenomicRegions() {
		return getReferenceGenomicRegions(new ArrayList<ImmutableReferenceGenomicRegion<D>>());
	}
	
	default Spliterator<ImmutableReferenceGenomicRegion<D>> iterateReferenceGenomicRegions(ReferenceSequence reference) {
		return immut(iterateMutableReferenceGenomicRegions(reference));
	}
	
	default Spliterator<? extends GenomicRegion> iterateGenomicRegions(ReferenceSequence reference) {
		return toRegion(iterateMutableReferenceGenomicRegions(reference));
	}
	
	default <C extends Collection<GenomicRegion>> C getRegionsIntersecting(ReferenceSequence reference, GenomicRegion region, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, region);
		split.forEachRemaining(mrgr->re.add(mrgr.getRegion()));
		return re;
	}
	
	default <C extends Map<GenomicRegion,D>> C getRegionsIntersecting(ReferenceSequence reference, GenomicRegion region, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, region);
		split.forEachRemaining(mrgr->re.put(mrgr.getRegion(),mrgr.getData()));
		return re;
	}
	
	default <C extends Collection<ImmutableReferenceGenomicRegion<D>>> C getReferenceRegionsIntersecting(ReferenceSequence reference, GenomicRegion region, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, region);
		split.forEachRemaining(mrgr->re.add(mrgr.toImmutable()));
		return re;
	}
	
	default <C extends Collection<GenomicRegion>> C getRegionsIntersecting(ReferenceSequence reference, int start, int end, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, start, end);
		split.forEachRemaining(mrgr->re.add(mrgr.getRegion()));
		return re;
	}
	
	default <C extends Map<GenomicRegion,D>> C getRegionsIntersecting(ReferenceSequence reference, int start, int end, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, start, end);
		split.forEachRemaining(mrgr->re.put(mrgr.getRegion(),mrgr.getData()));
		return re;
	}
	
	default <C extends Collection<ImmutableReferenceGenomicRegion<D>>> C getReferenceRegionsIntersecting(ReferenceSequence reference, int start, int end, C re) {
		Spliterator<? extends MutableReferenceGenomicRegion<D>> split = iterateIntersectingMutableReferenceGenomicRegions(reference, start, end);
		split.forEachRemaining(mrgr->re.add(mrgr.toImmutable()));
		return re;
	}
	
	
	// Concrete collection/map implementations
	default ArrayList<GenomicRegion> getRegionsIntersectingList(ReferenceSequence reference, GenomicRegion region) {
		return getRegionsIntersecting(reference, region, new ArrayList<GenomicRegion>());
	}
	
	default HashMap<GenomicRegion,D> getRegionsIntersectingMap(ReferenceSequence reference, GenomicRegion region) {
		return getRegionsIntersecting(reference, region, new HashMap<GenomicRegion,D>());
	}
	
	default ArrayList<ImmutableReferenceGenomicRegion<D>> getReferenceRegionsIntersecting(ReferenceSequence reference, GenomicRegion region) {
		return getReferenceRegionsIntersecting(reference, region, new ArrayList<ImmutableReferenceGenomicRegion<D>>());
	}
	
	default ArrayList<GenomicRegion> getRegionsIntersectingList(ReferenceSequence reference, int start, int end) {
		return getRegionsIntersecting(reference, start, end, new ArrayList<GenomicRegion>());
	}
	
	default HashMap<GenomicRegion,D> getRegionsIntersectingMap(ReferenceSequence reference, int start, int end) {
		return getRegionsIntersecting(reference, start, end, new HashMap<GenomicRegion,D>());
	}
	
	default ArrayList<ImmutableReferenceGenomicRegion<D>> getReferenceRegionsIntersecting(ReferenceSequence reference, int start, int end) {
		return getReferenceRegionsIntersecting(reference, start, end, new ArrayList<ImmutableReferenceGenomicRegion<D>>());
	}
	
	
	
	
	/**
	 * Changes get reflected into the storage, if supported
	 * @return
	 */
	default Collection<ImmutableReferenceGenomicRegion<D>> asCollection() {
		return new GenomicRegionStorageCollection<D>(this);
	}
	default D getRandomRecord() {
		Object[] re = new Object[1];
		iterateMutableReferenceGenomicRegions().tryAdvance(rgr->{re[0] = rgr.getData();});
		return (D) re[0];
	}

	default ReferenceGenomicRegion<D> getRandomEntry() {
		ReferenceGenomicRegion[] re = new ReferenceGenomicRegion[1];
		iterateMutableReferenceGenomicRegions().tryAdvance(rgr->{re[0] = rgr;});
		return re[0];
	}

	
	default <K,V> HashMap<K,V> index(Function<ReferenceGenomicRegion<D>,K> key, BiFunction<K,ReferenceGenomicRegion<D>,V> value) {
		HashMap<K, V> re = new HashMap<K, V>();
		iterateReferenceGenomicRegions().forEachRemaining(rgr->{
			K k = key.apply(rgr);
			V v = value.apply(k, rgr);
			re.put(k,v);
		});
		return re;
	}

	default String getName(){return "";}

	default String[] getMetaDataConditions() {
		if (getMetaData().hasProperty("conditions")) {
			int numCond = getMetaData().getEntry("conditions").asArray().length;
			String[] conditions = new String[numCond];
			for (int c=0; c<conditions.length; c++) {
				conditions[c] = getMetaData().getEntry("conditions").getEntry(c).getEntry("name").asString();
				if ("null".equals(conditions[c])) conditions[c] = c+"";
			}
			return conditions;
		}
		if (getMetaData().isNull() && getRandomRecord() instanceof HasConditions) {
			int numCond = ((HasConditions) getRandomRecord()).getNumConditions();
			String[] conditions = new String[numCond];
			for (int c=0; c<conditions.length; c++)
				conditions[c] = c+"";
			return conditions;
		}
		if (getMetaData().isNull() && getRandomRecord() instanceof NumericArray) {
			int numCond = ((NumericArray) getRandomRecord()).length();
			String[] conditions = new String[numCond];
			for (int c=0; c<conditions.length; c++)
				conditions[c] = c+"";
			return conditions;
		}
		return new String[0];
	}
	
	default double[] getMetaDataTotals() {
		if (getMetaData().hasProperty("conditions")) {
			int numCond = getMetaData().getEntry("conditions").asArray().length;
			double[] totals = new double[numCond];
			for (int c=0; c<totals.length; c++) 
				totals[c] = getMetaData().getEntry("conditions").getEntry(c).getEntry("total").asDouble();
			return totals;
		}
		return null;
	}
	
	default double[] getMetaDataTotals(String genome) {
		if (genome==null || genome.equals("")) return getMetaDataTotals();
		if (getMetaData().hasProperty("conditions")) {
			int numCond = getMetaData().getEntry("conditions").asArray().length;
			double[] totals = new double[numCond];
			for (int c=0; c<totals.length; c++) 
				totals[c] = getMetaData().getEntry("conditions").getEntry(c).getEntry("total_"+genome).asDouble();
			return totals;
		}
		return null;
	}

	public static <T> GenomicRegionStorage<T> load(String s) {
		try {
			Path p = Paths.get(s);
			return (GenomicRegionStorage<T>) WorkspaceItemLoaderExtensionPoint.getInstance().get(p).load(p);
		} catch (Exception e) {
			throw new RuntimeException("Could not load storage!",e);
		}
	}
	
}
