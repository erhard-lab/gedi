package gedi.slam;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.log1p;
import static java.lang.Math.max;

import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionMappingAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateFunctionPenaltyAdapter;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.MultiDirectionalSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleBounds;

import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import jdistlib.Beta;
import jdistlib.Binomial;

public class NumericalIntegrationProportion {

	
	private double alpha;
	private double beta;
	
	private double p;
	private double q;
	
	private double[] x;
	private double lower;
	private double upper;
	
	private boolean fitbeta;
	
	/**
	 * 
	 * @param alpha Prior
	 * @param beta Prior
	 * @param p Mismatch prob
	 * @param q Conversion prob
	 */
	public NumericalIntegrationProportion(double alpha, double beta, double p, double q, double lower, double upper, boolean fitbeta) {
		this.alpha = alpha;
		this.beta = beta;
		this.p = p;
		this.q = q;
		x = new double[101];
		for (int i=0; i<x.length; i++)
			x[i] = i/100.0;
		this.lower = lower;
		this.upper = upper;
		this.fitbeta = fitbeta;
	}
	
	
	public double[] infer(int[] d, int[] n, double[] m, Logger log) {
		double[] fs = new double[x.length];
		double[] fs2 = new double[x.length];
		double E = 0;
		
		for (int i=0; i<fs.length; i++) 
			fs[i] = loglikprior(x[i], d, n, m)-log(x.length-1);
		
		double[] p = fs.clone();
		
		fs[0]-=log(2);
		fs2[0] = fs[0];
		
		for (int i=1; i<fs.length; i++) { 
			fs2[i] =  lse(fs[i-1],fs[i]-log(2));
			fs[i] = lse(fs[i-1],fs[i]);
		}
		
		for (int i=0; i<fs.length; i++) 
			E+=x[i]*exp(p[i]-fs2[fs2.length-1]);
		
		// normalizing
		for (int i=0; i<fs.length; i++) 
			fs2[i] -= fs2[fs.length-1];
		
		for (int i=0; i<fs.length; i++)
			fs2[i] = exp(fs2[i]);
		
		int l = Arrays.binarySearch(fs2, lower);
		if (l<0) l=Math.max(0, -l-2);
		int u = Arrays.binarySearch(fs2, upper);
		if (u<0) u=-u-1;
		
		
		if (fitbeta) {
			// fit beta distribution
			MultivariateFunction fun = a->{
				double s = 0;
				for (int i=0; i<x.length; i++) {
					double ss = Beta.cumulative(x[i], a[0], a[1], true, false)-fs2[i];
					s+=ss*ss;
				}
				return s;
			};
			double[] beta = new SimplexOptimizer(1E-5,1E-5).optimize(
					GoalType.MINIMIZE,
					new NelderMeadSimplex(2),
//					new SimpleBounds(new double[]{1,1},new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY}),
					new ObjectiveFunction(fun),
					new InitialGuess(new double[]{1,1}),
					new MaxEval(1000)
					).getPointRef();
			
			if (beta[0]<=0 || beta[1]<=0)
				log.warning("Could not fit beta distribution!");
			
			return new double[] {x[l],Math.min(1, Math.max(0, Math.round(E*100)/100.0)),x[u],beta[0],beta[1]};
		}
		
		return new double[] {x[l],Math.min(1, Math.max(0, Math.round(E*100)/100.0)),x[u]};
	}
	
	
	private double loglikprior(double x, int[] d, int[] n, double[] m) {
		double s = 0;
		for (int i=0; i<d.length; i++) {
			double s1 = log(x)+Binomial.density(d[i], n[i], q, true);
			double s2 = log1p(-x)+Binomial.density(d[i], n[i], p, true);
			s+=lse(s1,s2)*m[i];
		}
		return s+Beta.density(x, alpha, beta, true);
	}
	
	private static final double lse(double u, double v) {
	    double m = max(u,v);
	    return log(exp(u-m)+exp(v-m))+m;
	}
	
	private static final double lse(double[] x) {
		double m = ArrayUtils.max(x);
		double s = 0;
		for (double xx : x)
			s+=exp(xx-m);
		return log(s)+m;
	}
	
}
