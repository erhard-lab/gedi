package gedi.util.datastructure.collections;

import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.orm.OrmSerializer;

import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;

/**
 * Not thread-safe!
 * @author erhard
 *
 * @param <T>
 */
public class SortingCollection<T> implements Collection<T>, Closeable {

	private int memoryCapacity;
	private ArrayList<T> mem;
	private PageFileWriter buf = null;
	private ArrayList<Long> offsets = new ArrayList<Long>(); 
	private Class<T> cls;
	private Comparator<? super T> comp;
	private int size;
	
	private OrmSerializer orm = new OrmSerializer(true,false);
	
	public SortingCollection(Class<T> cls, Comparator<? super T> comp, int memoryCapacity) {
		this.cls = cls;
		this.comp = comp;
		this.memoryCapacity = memoryCapacity;
		mem = new ArrayList<T>(memoryCapacity);
	}
	
	
	public boolean add(T e) {
		size++;
		if (mem.size()>=memoryCapacity) 
			spillToDisk();
		mem.add(e);
		return true;
	}

	
	public ExtendedIterator<T> iterator() {
		if (buf==null) {
			Collections.sort(mem, comp);
			return EI.wrap(mem);
		}
		
		spillToDisk();
		try {
			PageFile f = buf.read(false);
			ExtendedIterator<T>[] iter = new ExtendedIterator[offsets.size()];
			for (int i=0; i<iter.length; i++) {
				iter[i] = (ExtendedIterator<T>) f.view(i==0?0:offsets.get(i-1), offsets.get(i)).iterator(p->{
					try {
						T re = (T)orm.deserialize(p);
						orm.clearObjectCache();
						return re;
					} catch (Exception e) {
						throw new RuntimeException("Could not deserialize object!",e);
					}
				}
				);
			}
			return EI.merge(comp, iter);
		
		} catch (IOException e) {
			throw new RuntimeException("Cannot iterate temp file!",e);
		}
	}
	
	
	@Override
	public void close() throws IOException {
		if (buf!=null) {
			buf.close();
			new File(buf.getPath()).delete();
			buf = null;
		}
		mem.clear();
	}

	private void spillToDisk() {
		Collections.sort(mem, comp);
		if (buf==null)
			try {
				buf = new PageFileWriter(Files.createTempFile("sortingcollection", ".tmp").toString());
				new File(buf.getPath()).deleteOnExit();
			} catch (IOException e) {
				throw new RuntimeException("Cannot create temp file!",e);
			}
		
				
		for (T e : mem)
			try {
				orm.serialize(buf, e);
				orm.clearObjectCache();
			} catch (IOException e1) {
				throw new RuntimeException("Cannot write entry !",e1);
			}
		
		offsets.add(buf.position());
		
		mem.clear();
	}


	@Override
	public int size() {
		return size;
	}


	@Override
	public boolean isEmpty() {
		return size==0;
	}


	@Override
	public boolean contains(Object o) {
		throw new UnsupportedOperationException();
	}


	@Override
	public Object[] toArray() {
		return iterator().toArray(cls);
	}


	@Override
	public <E> E[] toArray(E[] a) {
		return (E[]) iterator().toArray((T[]) a);
	}


	@Override
	public boolean remove(Object o) {
		throw new UnsupportedOperationException();
	}


	@Override
	public boolean containsAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}


	@Override
	public boolean addAll(Collection<? extends T> c) {
		for (T t : c) add(t);
		return true;
	}


	@Override
	public boolean removeAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}


	@Override
	public boolean retainAll(Collection<?> c) {
		throw new UnsupportedOperationException();
	}


	@Override
	public void clear() {
		try {
			close();
		} catch (IOException e) {
			throw new RuntimeException("Could not clear!",e);
		}
	}
	
	
}
