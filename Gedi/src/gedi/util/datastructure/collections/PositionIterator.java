package gedi.util.datastructure.collections;

import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.functions.ExtendedIterator;

public interface PositionIterator<T> extends IntIterator {

	T getData();
	
	default T nextData() {
		nextInt();
		return getData();
	}
	
	default ExtendedIterator<T> data() {
		return new ExtendedIterator<T>() {

			@Override
			public boolean hasNext() {
				return PositionIterator.this.hasNext();
			}

			@Override
			public T next() {
				return PositionIterator.this.nextData();
			}
			
		};
	}
	
}
