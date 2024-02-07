package gedi.util.datastructure.collections.intcollections;

import gedi.util.FunctorUtils;
import gedi.util.FunctorUtils.DemultiplexIntIterator;
import gedi.util.FunctorUtils.DemultiplexIterator;
import gedi.util.FunctorUtils.FilteredIntIterator;
import gedi.util.FunctorUtils.MappedIterator;
import gedi.util.FunctorUtils.PeekIntIterator;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.ExtendedIterator;

import java.util.Iterator;
import java.util.Objects;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.IntBinaryOperator;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntPredicate;
import java.util.function.IntUnaryOperator;

public interface IntIterator extends ExtendedIterator<Integer> {
	
	
	default IntIterator filterInt(IntPredicate predicate) {
		return new FilteredIntIterator(this,predicate);
	}
	
	default <T> ExtendedIterator<T> mapInt(IntFunction<T> mapper) {
		return new FunctorUtils.MappedIntIterator<>(this,mapper);
	}
	
	
	default IntIterator mapIntToInt(IntUnaryOperator mapper) {
		return new FunctorUtils.MappedIntToIntIterator(this,mapper);
	}
	
	default <O> DemultiplexIntIterator<O> unfoldInt(IntFunction<Iterator<O>> demulti) {
		return FunctorUtils.demultiplexIterator(this, demulti);
	}
	
	
	default PeekIntIterator peekInt() {
		return new PeekIntIterator(this);
	}

    default void forEachRemainingInt(IntConsumer action) {
        Objects.requireNonNull(action);
        while (hasNext())
            action.accept(nextInt());
    }
    
	default Integer next() {
		return nextInt();
	}
	
	public int nextInt();
	

	default NumericArray toNumericArray() {
		return NumericArray.wrap(toIntArray());
	}

	
	public static final IntIterator empty = new EmptyIntIterator();
	public static IntIterator singleton(int value) { return new SingletonIntIterator(value);}
	
	
	public static class EmptyIntIterator implements IntIterator {
		@Override
		public int nextInt() {
			return -1;
		}

		@Override
		public boolean hasNext() {
			return false;
		}

		@Override
		public Integer next() {
			return null;
		}

		@Override
		public void remove() {
		}
		
	}
	
	public static class SingletonIntIterator implements IntIterator {
		private int v;
		private boolean fresh = true;
		
		public SingletonIntIterator(int v) {
			this.v = v;
		}

		@Override
		public boolean hasNext() {
			return fresh;
		}

		@Override
		public Integer next() {
			fresh = false;
			return v;
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}

		@Override
		public int nextInt() {
			fresh = false;
			return v;
		}

	}
	
	public static class ArrayIterator implements IntIterator {
		private int next;
		private int[] a;
		private int end;
		
		public ArrayIterator(int[] a) {
			this(a,0,a.length);
		}
		public ArrayIterator(int[] a,int start, int end) {
			this.a = a;
			this.end = end;
			this.next = start;
		}

		@Override
		public boolean hasNext() {
			return next<end;
		}

		@Override
		public Integer next() {
			return a[next++];
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}

		@Override
		public int nextInt() {
			return a[next++];
		}

	}
	
	
	default int sum() {
		int  s = 0;
		while(hasNext())
			s+=nextInt();
		return s;
	}
	
	default int max() {
		return stat(Integer.MIN_VALUE,Math::max);
	}
	
	default int min() {
		return stat(Integer.MAX_VALUE,Math::min);
	}
	
	default int absmax() {
		return stat(0,(a,b)->Math.abs(a)>Math.abs(b)?a:b);
	}
	
	default int absmin() {
		return stat(Integer.MAX_VALUE,(a,b)->Math.abs(a)<Math.abs(b)?a:b);
	}

	
	default int stat(int init, IntBinaryOperator op) {
		int s = init;
		while(hasNext())
			s=op.applyAsInt(s, nextInt());
		return s;
	}
	
}
