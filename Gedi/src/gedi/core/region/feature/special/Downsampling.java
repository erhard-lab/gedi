package gedi.core.region.feature.special;

import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.functions.NumericArrayFunction;

public enum Downsampling {

	Max {
		@Override
		public NumericArray downsample(NumericArray a) {
			double max = a.evaluate(NumericArrayFunction.Max);
			if (max<=1) return a;
			for (int i=0; i<a.length(); i++)
				a.setDouble(i,a.getDouble(i)/max);
			return a;
		}
	},

	Logsc {
		@Override
		public NumericArray downsample(NumericArray a) {
			double max = a.evaluate(NumericArrayFunction.Max);
			if (max<=1) return a;
			double fac = Math.log(max)/Math.log(2)/max;
			for (int i=0; i<a.length(); i++)
				a.setDouble(i,a.getDouble(i)*fac);
			return a;
		}
	},
	
	No{
		@Override
		public NumericArray downsample(NumericArray a) {
			return a;
		}
	},
	
	Digital {
		@Override
		public NumericArray downsample(NumericArray a) {
			for (int i=0; i<a.length(); i++)
				if (a.getDouble(i)>1)
					a.setDouble(i,1);
			return a;
		}
	};
	
	public abstract NumericArray downsample(NumericArray a);
	
	
}
