package gedi.util.math.stat.testing;

import gedi.util.math.stat.DoubleRanking;

public class MultipleTestingCorrection {

	
	/**
	 * Inplace!!
	 * @param p
	 * @return
	 */
	public static double[] benjaminiHochberg(double[] p) {
		double n = p.length;
		DoubleRanking r = new DoubleRanking(p);
		r.sort(false);
		double min = 1;
		for (int i=0; i<p.length; i++) {
			min = p[i] = Math.min(min, n/(n-i) * p[i]);
		}
		r.restore();
		return p;
	}
	
}
