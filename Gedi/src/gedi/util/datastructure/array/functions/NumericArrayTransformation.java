package gedi.util.datastructure.array.functions;

import java.util.Random;
import java.util.function.Consumer;
import java.util.function.DoubleUnaryOperator;
import java.util.function.UnaryOperator;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.math.stat.DoubleRanking;
import gedi.util.math.stat.RandomNumbers;

/**
 * In place!!!
 * @author erhard
 *
 */
public interface NumericArrayTransformation extends UnaryOperator<NumericArray> {

	
	public static NumericArrayTransformation Rank = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			DoubleRanking dr = new DoubleRanking(t.toDoubleArray());
			dr.sort(true);
			for (int i=0; i<dr.size(); i++)
				t.setInt(dr.getOriginalIndex(i), i);
			return t;
		}
	};
	
	public static NumericArrayTransformation RelativeRank = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			DoubleRanking dr = new DoubleRanking(t.toDoubleArray());
			dr.sort(true);
			for (int i=0; i<dr.size(); i++)
				t.setDouble(dr.getOriginalIndex(i), i/(double)(dr.size()-1));
			return t;
		}
	};
	
	public static NumericArrayTransformation InverseRank = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			DoubleRanking dr = new DoubleRanking(t.toDoubleArray());
			dr.sort(false);
			for (int i=0; i<dr.size(); i++)
				t.setInt(dr.getOriginalIndex(i), i);
			return t;
		}
	};
	
	public static NumericArrayTransformation InverseRelativeRank = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			DoubleRanking dr = new DoubleRanking(t.toDoubleArray());
			dr.sort(false);
			for (int i=0; i<dr.size(); i++)
				t.setDouble(dr.getOriginalIndex(i), i/(double)(dr.size()-1));
			return t;
		}
	};
	
	public static NumericArrayTransformation CumSum = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			for (int i=1; i<t.length(); i++)
				t.setDouble(i,t.getDouble(i-1)+t.getDouble(i));
			return t;
		}
	};
	
	public static NumericArrayTransformation InverseCumSum = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			for (int i=t.length()-2; i>=0; i--)
				t.setDouble(i,t.getDouble(i+1)+t.getDouble(i));
			return t;
		}
	};
	
	public static NumericArrayTransformation Unitize = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			double max = t.evaluate(NumericArrayFunction.Max);
			for (int i=0; i<t.length(); i++)
				t.setDouble(i,t.getDouble(i)/max);
			return t;
		}
	};
	
	public static NumericArrayTransformation SubtractFrom(double from) { 
		return (t)-> {
			for (int i=0; i<t.length(); i++) {
				t.add(i, -from);
				t.mult(i, -1);
			}
			return t;
		};
	}
	

	
	public static NumericArrayTransformation add(double value) { 
		return (t)-> {
			for (int i=0; i<t.length(); i++) {
				t.add(i, value);
			}
			return t;
		};
	}
	
	public static NumericArrayTransformation constant(double value) { 
		return (t)-> {
			for (int i=0; i<t.length(); i++) {
				t.setDouble(i, value);
			}
			return t;
		};
	}

	
	public static NumericArrayTransformation mult(double value) { 
		return (t)-> {
			for (int i=0; i<t.length(); i++) {
				t.mult(i, value);
			}
			return t;
		};
	}
	
	public static NumericArrayTransformation reciprocal() { 
		return function(n->1/n);
	}
	
	public static NumericArrayTransformation function(DoubleUnaryOperator op) { 
		return (t)-> {
			for (int i=0; i<t.length(); i++) {
				t.setDouble(i, op.applyAsDouble(t.getDouble(i)));
			}
			return t;
		};
	}

	public static NumericArrayTransformation Invert = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			NumericArray buff = NumericArray.createMemory(1, t.getType());
			for (int i=0; i<t.length()/2; i++) {
				buff.copy(t, i, 0);
				t.copy(t.length()-1-i, i);
				t.copy(buff,0,t.length()-1-i);
			}
			return t;
		}
	};
	
	public static NumericArrayTransformation Shuffle = new NumericArrayTransformation() {
		@Override
		public NumericArray apply(NumericArray t) {
			RandomNumbers rnd = RandomNumbers.getGlobal();
			NumericArray buff = NumericArray.createMemory(1, t.getType());
			for (int k=0; k<t.length()-1; k++) {
				int m =rnd.getUnif(0,t.length()-k)+k;
				buff.copy(t, k, 0);
				t.copy(m, k);
				t.copy(buff,0,m);
			}
			return t;
		}
	};

	
	
	
}
