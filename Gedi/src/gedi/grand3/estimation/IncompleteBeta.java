package gedi.grand3.estimation;

import org.apache.commons.math3.exception.ConvergenceException;
import org.apache.commons.math3.exception.MaxCountExceededException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Precision;

public class IncompleteBeta {

	
	private static double b(int n, double x, double a, double b) {
		double ret;
		double m;
		if (n % 2 == 0) { // even
			m = n / 2.0;
			ret = (m * (b - m) * x) / ((a + (2 * m) - 1) * (a + (2 * m)));
		} else {
			m = (n - 1.0) / 2.0;
			ret = -((a + m) * (a + b + m) * x) / ((a + (2 * m)) * (a + (2 * m) + 1.0));
		}
		return ret;
	}
	
	public static double libeta(double x, final double a, final double b) {
//		if (x > (a + 1) / (2 + b + a) && 1 - x <= (b + 1) / (2 + b + a)) 
//			throw new RuntimeException(x+" "+a+" "+b);
//		return (a * FastMath.log(x)) + (b * FastMath.log1p(-x)) - FastMath.log(a) 
//				- Math.log(cfrac(x, a,b));
		return FastMath.log(hyperg_2F1_series(a+b, 1, a+1, x)) + (a * FastMath.log(x)) + (b * FastMath.log1p(-x)) - FastMath.log(a);
	}
	
	/**
	 * This is directly translated from GSL
	 */
	private static final double GSL_DBL_EPSILON=2.2204460492503131e-16;
	public static double hyperg_2F1_series(double a, double b, double c,
            double x)
		{
		double sum_pos = 1.0;
		double sum_neg = 0.0;
		double del_pos = 1.0;
		double del_neg = 0.0;
		double del = 1.0;
		double del_prev;
		double k = 0.0;
		int i = 0;
		
		if(Math.abs(c) < GSL_DBL_EPSILON) {
		return Double.NaN;
		}
		
		do {
		if(++i > 30000) {
			return Double.NaN;
		}
		del_prev = del;
		del *= (a+k)*(b+k) * x / ((c+k) * (k+1.0));  /* Gauss series */
		
		if(del > 0.0) {
		del_pos  =  del;
		sum_pos +=  del;
		}
		else if(del == 0.0) {
		/* Exact termination (a or b was a negative integer).
		 */
		del_pos = 0.0;
		del_neg = 0.0;
		break;
		}
		else {
		del_neg  = -del;
		sum_neg -=  del;
		}
		
		/*
		* This stopping criteria is taken from the thesis
		* "Computation of Hypergeometic Functions" by J. Pearson, pg. 31
		* (http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf)
		* and fixes bug #45926
		*/
		if (Math.abs(del_prev / (sum_pos - sum_neg)) < GSL_DBL_EPSILON &&
				Math.abs(del / (sum_pos - sum_neg)) < GSL_DBL_EPSILON)
		break;
		
		k += 1.0;
		} while(Math.abs((del_pos + del_neg)/(sum_pos-sum_neg)) > GSL_DBL_EPSILON);
		
		return sum_pos - sum_neg;
	}
	
	private static double cfrac(double x, double pa, double pb) {
		final double small = 1e-50;
        double hPrev = 1;
        int n = 1;
        double dPrev = 0.0;
        double cPrev = hPrev;
        double hN = hPrev;

        while (n < Integer.MAX_VALUE) {
            final double a = 1;
            final double b = b(n, x,pa,pb);

            double dN = a + b * dPrev;
            if (Precision.equals(dN, 0.0, small)) {
                dN = small;
            }
            double cN = a + b / cPrev;
            if (Precision.equals(cN, 0.0, small)) {
                cN = small;
            }

            dN = 1 / dN;
            final double deltaN = cN * dN;
            hN = hPrev * deltaN;

            if (Double.isInfinite(hN)) {
                throw new ConvergenceException(LocalizedFormats.CONTINUED_FRACTION_INFINITY_DIVERGENCE,
                                               x);
            }
            if (Double.isNaN(hN)) {
                throw new ConvergenceException(LocalizedFormats.CONTINUED_FRACTION_NAN_DIVERGENCE,
                                               x);
            }

            if (FastMath.abs(deltaN - 1.0) < 1E-14) {
                break;
            }

            dPrev = dN;
            cPrev = cN;
            hPrev = hN;
            n++;
        }

        if (n >= Integer.MAX_VALUE) {
            throw new MaxCountExceededException(LocalizedFormats.NON_CONVERGENT_CONTINUED_FRACTION,
            		Integer.MAX_VALUE, x);
        }

        return hN;
	}
}
