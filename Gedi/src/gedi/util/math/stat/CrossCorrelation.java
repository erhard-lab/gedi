package gedi.util.math.stat;

import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import gedi.util.ArrayUtils;

public class CrossCorrelation {

	
	public double[] padAndcompute(double[] a, double[] b) {
		a = ArrayUtils.padToPowerOfTwo(a);
		b = ArrayUtils.padToPowerOfTwo(b);
		if (a.length<b.length) a = ArrayUtils.redimPreserve(a, b.length);
		if (b.length<a.length) b = ArrayUtils.redimPreserve(b, a.length);
		return compute(a,b);
	}
	public double[] compute(double[] a, double[] b) {
		
		if (a.length!=b.length) {
			throw new RuntimeException("Must be of same length!");
		}
		
				
		double[][] ah = {
				a,
				new double[a.length]
		};
		double[][] bh = {
				b,
				new double[b.length]
		};
		FastFourierTransformer.transformInPlace(ah, DftNormalization.STANDARD, TransformType.FORWARD);
		FastFourierTransformer.transformInPlace(bh, DftNormalization.STANDARD, TransformType.FORWARD);
		for (int p=0; p<a.length; p++) {
			double real = ah[0][p]*bh[0][p]+ah[1][p]*bh[1][p];
			double imag = ah[0][p]*bh[1][p]-bh[0][p]*ah[1][p];
			ah[0][p] = real;
			ah[1][p] = imag;
		}
		FastFourierTransformer.transformInPlace(ah, DftNormalization.STANDARD, TransformType.INVERSE);
		return ah[0];
	}
	
}
