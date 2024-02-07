package gedi.util.datastructure.collections.doublecollections;

import java.util.function.DoublePredicate;

import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;

import gedi.util.FunctorUtils.FilteredDoubleIterator;
import gedi.util.MathUtils;
import gedi.util.functions.ExtendedIterator;

public interface DoubleIterator extends ExtendedIterator<Double> {

	
	default DoubleIterator filterDouble(DoublePredicate predicate) {
		return new FilteredDoubleIterator(this,predicate);
	}
	
	default Double next() {
		return nextDouble();
	}
	
	
	public double nextDouble();

	public static class EmptyDoubleIterator implements DoubleIterator {
		@Override
		public double nextDouble() {
			return -1;
		}

		@Override
		public boolean hasNext() {
			return false;
		}

		@Override
		public Double next() {
			return null;
		}

		@Override
		public void remove() {
		}
		
	}
	
	public static class ArrayIterator implements DoubleIterator {
		private int next;
		private double[] a;
		private int end;
		
		public ArrayIterator(double[] a) {
			this(a,0,a.length);
		}
		public ArrayIterator(double[] a,int start, int end) {
			this.a = a;
			this.end = end;
			this.next = start;
		}

		@Override
		public boolean hasNext() {
			return next<end;
		}

		@Override
		public Double next() {
			return a[next++];
		}

		@Override
		public void remove() {
			throw new UnsupportedOperationException();
		}

		@Override
		public double nextDouble() {
			return a[next++];
		}

	}
	
	default double evaluate(UnivariateStatistic stat) {
		return stat.evaluate(toDoubleArray());
	}
	
	default double sum() {
		double s = 0;
		while(hasNext())
			s+=nextDouble();
		return s;
	}
	
	default double mean() {
		double s = 0;
		int n = 0;
		while(hasNext()) {
			s+=nextDouble();
			n++;
		}
		return s/n;
	}
	
	default int argmin() {
		int index = 0;
		int minIndex = 0;
		double s = Double.POSITIVE_INFINITY;
		while(hasNext()) {
			double n = nextDouble();
			if (n<s) {
				s = n;
				minIndex = index;
			}
			index++;
		}
		return minIndex;
	}
	
	default int argmax() {
		int index = 0;
		int maxIndex = 0;
		double s = Double.NEGATIVE_INFINITY;
		while(hasNext()) {
			double n = nextDouble();
			if (n>s) {
				s = n;
				maxIndex = index;
			}
			index++;
		}
		return maxIndex;
	}
	
	default double min() {
		double s = Double.POSITIVE_INFINITY;
		while(hasNext())
			s=Math.min(s, nextDouble());
		return s;
	}
	
	default double max() {
		double s = Double.NEGATIVE_INFINITY;
		while(hasNext())
			s=Math.max(s, nextDouble());
		return s;
	}
	
	default double saveMax(double allUnsafeValue) {
		double s = Double.NEGATIVE_INFINITY;
		while(hasNext())
			s=MathUtils.saveMax(s, nextDouble());
		if (Double.isInfinite(s))return allUnsafeValue;
		return s;
	}
	
	default double saveMin(double allUnsafeValue) {
		double s = Double.POSITIVE_INFINITY;
		while(hasNext())
			s=MathUtils.saveMin(s, nextDouble());
		if (Double.isInfinite(s))return allUnsafeValue;
		return s;
	}
	
	default double[] saveMinMax(double allUnsafeValueMin, double allUnsafeValueMax) {
		double min = Double.POSITIVE_INFINITY;
		double max = Double.NEGATIVE_INFINITY;
		while(hasNext()) {
			double v = nextDouble();
			min=MathUtils.saveMin(min, v);
			max=MathUtils.saveMax(max, v);
		}
		return new double[] {
				Double.isInfinite(min)?allUnsafeValueMin:min,
				Double.isInfinite(max)?allUnsafeValueMax:max
		};
	}
	
}
