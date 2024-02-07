package gedi.core.data.reads;

import java.util.HashMap;
import java.util.function.DoubleToIntFunction;
import java.util.function.DoubleUnaryOperator;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.datastructure.array.sparse.AutoSparseDenseIntArrayCollector;

public interface ReadCountMode {

	
	public static final ReadCountMode All = new ReadCountMode() {
		
		public double computeCount(int count, int multiplicity, double weight) {
			return count;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			return count;
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return count;
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count);
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count);
		}
		
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count);
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setInt(destIndex,count); 
		}
		
		public Class<? extends Number> getType() {
			return Integer.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Integer;
		}
		public String name() { return "All"; }
		@Override
		public String toString() {
			return name();
		}
		
	};
	
	public static final ReadCountMode Weight = new ReadCountMode() {
		
		public double computeCount(int count, int multiplicity, double weight) {
			return count*weight;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			throw new IllegalArgumentException("Mode weight produces real counts!");
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return (int)(count*weight);
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count*weight); 
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count*weight); 
		}
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			throw new IllegalArgumentException("Mode weight produces real counts!");
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setDouble(destIndex,count*weight); 
		}
		
		public Class<? extends Number> getType() {
			return Double.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Double;
		}
		
		public String name() { return "Weight"; }
		@Override
		public String toString() {
			return name();
		}
		
	};
	
	public static final ReadCountMode Divide = new ReadCountMode() {
		public double computeCount(int count, int multiplicity, double weight) {
			return (double)count/multiplicity;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			throw new IllegalArgumentException("Mode divide produces real counts!");
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return (int)((double)count/multiplicity);
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,(double)count/multiplicity);
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,(double)count/multiplicity);
		}
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			throw new IllegalArgumentException("Mode weight produces real counts!");
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setDouble(destIndex,(double)count/multiplicity);
		}
		
		public Class<? extends Number> getType() {
			return Double.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Double;
		}
		public String name() { return "Divide"; }
		@Override
		public String toString() {
			return name();
		}
		
	};

	public static final ReadCountMode Unique = new ReadCountMode() {
		public double computeCount(int count, int multiplicity, double weight) {
			return multiplicity<=1?count:0;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			return multiplicity<=1?count:0;
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return multiplicity<=1?count:0;
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,multiplicity<=1?count:0);
		}
		
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,multiplicity<=1?count:0);
		}
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,multiplicity<=1?count:0);
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setInt(destIndex,multiplicity<=1?count:0);
		}
		
		public Class<? extends Number> getType() {
			return Integer.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Integer;
		}
		
		public String name() { return "Unique"; }
		@Override
		public String toString() {
			return name();
		}
		
	};

	public static final ReadCountMode CollapseAll = new ReadCountMode() {
		public double computeCount(int count, int multiplicity, double weight) {
			return count>0?1:0;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			return count>0?1:0;
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return count>0?1:0;
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0?1:0);
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0?1:0);
		}
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0?1:0);
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setInt(destIndex,count>0?1:0);
		}
		
		public Class<? extends Number> getType() {
			return Integer.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Integer;
		}
		
		public String name() { return "Collapse"; }
		@Override
		public String toString() {
			return name();
		}
		
	};

	public static final ReadCountMode CollapseUnique = new ReadCountMode() {
		public double computeCount(int count, int multiplicity, double weight) {
			return count>0&&multiplicity<=1?1:0;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			return count>0&&multiplicity<=1?1:0;
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return count>0&&multiplicity<=1?1:0;
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0&&multiplicity<=1?1:0);
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0&&multiplicity<=1?1:0);
		}
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
			dest.add(destIndex,count>0&&multiplicity<=1?1:0);
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
			dest.setInt(destIndex,count>0&&multiplicity<=1?1:0);
		}
		
		public Class<? extends Number> getType() {
			return Integer.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Integer;
		}
		
		public String name() { return "CollapseUnique"; }
		@Override
		public String toString() {
			return name();
		}
		
	};

	
	public default ReadCountMode transformCounts(double factor) {
		if (factor==0) return ReadCountMode.No;
		if (factor==1) return this;
		return transformCounts(x->x*factor);
	}
	public default ReadCountMode transformCounts(DoubleUnaryOperator transformer) {
		return new ReadCountMode() {
			public double computeCount(int count, int multiplicity, double weight) {
				return transformer.applyAsDouble(ReadCountMode.this.computeCount(count,multiplicity,weight));
			}
			
			public int computeCountInt(int count, int multiplicity, double weight) {
				throw new IllegalArgumentException("Transformed mode produces real counts!");
			}
			
			public int computeCountFloor(int count, int multiplicity, double weight) {
				throw new IllegalArgumentException("Transformed mode produces real counts!");
			}

			public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
				if (dest.isIntegral()) throw new IllegalArgumentException("Transformed mode produces real counts!");
				dest.add(destIndex, computeCount(count, multiplicity, weight));
			}
			
			public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
				dest.add(destIndex, computeCount(count, multiplicity, weight));
			}
			
			public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
				throw new IllegalArgumentException("Transformed mode produces real counts!");
			}
			
			public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
				if (dest.isIntegral()) throw new IllegalArgumentException("Transformed mode produces real counts!");
				dest.setDouble(destIndex, computeCount(count, multiplicity, weight));
			}
			
			public Class<? extends Number> getType() {
				return Double.TYPE;
			}
			
			public NumericArrayType getNumericArrayType() {
				return NumericArrayType.Double;
			}
			
			public String name() { return "Transformed"; }
			@Override
			public String toString() {
				return name();
			}
		};
	}
	
	
	
	public static final ReadCountMode No = new ReadCountMode() {
		
		public double computeCount(int count, int multiplicity, double weight) {
			return 0;
		}
		
		public int computeCountInt(int count, int multiplicity, double weight) {
			return 0;
		}
		
		public int computeCountFloor(int count, int multiplicity, double weight) {
			return 0;
		}

		public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
		}
		
		public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
		}
		
		public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
		}
		
		public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
		}
		
		public Class<? extends Number> getType() {
			return Integer.TYPE;
		}
		
		public NumericArrayType getNumericArrayType() {
			return NumericArrayType.Integer;
		}
		public String name() { return "No"; }
		@Override
		public String toString() {
			return name();
		}
		
	};
	public default ReadCountMode transformCountsInt(DoubleToIntFunction transformer) {
		return new ReadCountMode() {
			public double computeCount(int count, int multiplicity, double weight) {
				return transformer.applyAsInt(computeCount(count,multiplicity,weight));
			}
			
			public int computeCountInt(int count, int multiplicity, double weight) {
				return transformer.applyAsInt(computeCount(count,multiplicity,weight));
			}
			
			public int computeCountFloor(int count, int multiplicity, double weight) {
				return transformer.applyAsInt(computeCount(count,multiplicity,weight));
			}

			public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
				dest.add(destIndex, computeCountInt(count, multiplicity, weight));
			}
			
			public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
				dest.add(destIndex, computeCountInt(count, multiplicity, weight));
			}
			
			public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight) {
				dest.add(destIndex, computeCountInt(count, multiplicity, weight));
			}
			
			public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight) {
				dest.setDouble(destIndex, computeCountInt(count, multiplicity, weight));
			}
			
			public Class<? extends Number> getType() {
				return Integer.TYPE;
			}
			
			public NumericArrayType getNumericArrayType() {
				return NumericArrayType.Integer;
			}
			
			public String name() { return "Transformed"; }
			@Override
			public String toString() {
				return name();
			}
		};
	}
	
	
	
	public double computeCount(int count, int multiplicity, double weight) ;
	
	public int computeCountInt(int count, int multiplicity, double weight);
	
	public int computeCountFloor(int count, int multiplicity, double weight);

	public void addCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight);
	
	public void addCount(AutoSparseDenseDoubleArrayCollector dest, int destIndex, int count, int multiplicity, double weight);
	
	public void addCountInt(AutoSparseDenseIntArrayCollector dest, int destIndex, int count, int multiplicity, double weight);
	
	public void getCount(NumericArray dest, int destIndex, int count, int multiplicity, double weight);
	
	public Class<? extends Number> getType();
	
	public NumericArrayType getNumericArrayType();

	public String name();
	
	
	public static final ReadCountMode[] values = {
		All,Weight,Divide,Unique,CollapseAll,CollapseUnique
	};
	
	public static final HashMap<String,ReadCountMode> valueOf = ArrayUtils.index(values,v->v.name());
		
	public static ReadCountMode[] values() {
		return values;
	}

	public static ReadCountMode valueOf(String name) {
		return valueOf.get(name);
	}

	

	
	
}
