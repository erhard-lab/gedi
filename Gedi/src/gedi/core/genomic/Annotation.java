package gedi.core.genomic;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.function.Predicate;

import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;

public class Annotation<T> {

	private String id;
	
	private LinkedList<GenomicRegionStorage<T>> storages = new LinkedList<GenomicRegionStorage<T>>();;
	private MemoryIntervalTreeStorage<T> mem;
	private Predicate<ImmutableReferenceGenomicRegion<T>> pred;
	
	public Annotation(String id) {
		this.id = id;
	}
	
	public Annotation<T> set(GenomicRegionStorage<T> storage) {
		if (this.storages.size()>0) throw new IllegalArgumentException("Already set!");
		this.storages.add(storage);
		return this;
	}
	
	public void filter(Predicate<ImmutableReferenceGenomicRegion<T>> pred) {
		if (mem!=null) throw new RuntimeException("Already loaded!");
		this.pred = pred;
	}
	
	
	public void merge(Annotation<T> other) {
		if (mem!=null) {
			if (other.mem!=null)
				addToMem(other.mem);
			else {
				Iterator<GenomicRegionStorage<T>> it = other.storages.iterator();
				while (it.hasNext()) {
					MemoryIntervalTreeStorage<T> o = it.next().toMemory();
					addToMem(o);
				}
			}
		}
		else {
			if (other.mem!=null)
				storages.add(other.mem);
			else
				storages.addAll(other.storages);
		}
		
	}
	
	public MemoryIntervalTreeStorage<T> get() {
		if (mem==null) {
			Iterator<GenomicRegionStorage<T>> it = storages.iterator();
			while (it.hasNext()) {
				MemoryIntervalTreeStorage<T> o = it.next().toMemory();
				if (mem==null) mem = new MemoryIntervalTreeStorage<T>(o.getType());
				addToMem(o);
			}
		}
		return mem;
	}


	private void addToMem(MemoryIntervalTreeStorage<T> o) {
		long before = mem.size();
		long add = o.size();
		if (pred!=null)
			mem.fill(o.ei().filter(pred));
		else
			mem.fill(o);
		if (mem.size()<before+add && pred==null)
			throw new RuntimeException("Removed multi entries!");		
	}

	public String getId() {
		return id;
	}

}
