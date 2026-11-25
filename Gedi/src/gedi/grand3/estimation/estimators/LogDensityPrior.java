package gedi.grand3.estimation.estimators;

import jdistlib.Normal;

public enum LogDensityPrior {
		Epanechnikov {
			public double logdens(double x, double mean, double rad) {
				double nd = (x-mean)/rad;
				nd = nd*nd;
				if (nd>=1) return Double.NEGATIVE_INFINITY;
//				return Math.log(0.75*(1-nd)/rad);
				return Math.log1p(-nd);
			}
		},
		Gaussian {
			public double logdens(double x, double mean, double rad) {
				return Normal.density(x, mean, rad/2, true);
			}
		};
		
		public abstract double logdens(double x, double mean, double rad);
	}