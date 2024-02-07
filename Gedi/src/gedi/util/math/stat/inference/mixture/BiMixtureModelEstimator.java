package gedi.util.math.stat.inference.mixture;


import java.util.Arrays;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.UnivariateSolverUtils;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;

import jdistlib.Beta;
import net.jafama.FastMath;

public class BiMixtureModelEstimator {

	private BiMixtureModelData[] data;
//	private double priorAlpha = 1;
//	private double priorBeta = 1;

	

	public BiMixtureModelEstimator(BiMixtureModelData[] data) {
		this.data = data;
	}
	
//	public BiMixtureModel(BiMixtureModelData[] data, double priorAlpha, double priorBeta) {
//		this.data = data;
//		this.priorAlpha = priorAlpha;
//		this.priorBeta = priorBeta;
//	}
	
	/**
	 * if ll0>>ll1 in the data, the map will be close to 1 
	 * @param ci
	 * @param betaApprox
	 * @return
	 */
	public BiMixtureModelResult estimate(double ci, boolean betaApprox) {
		
		double lower = (1-ci)/2;
		double upper = 1-lower;
		
		double max = computeMAP();
		double ymax = loglik(max);
		
//		double left = bisect_inc(0,max,ymax+log(1E-3),1000,1E-6);
//		double right = bisect_dec(max,1,ymax+log(1E-3),1000,1E-6);
		
		double left = bisect(0,max,ymax+FastMath.log(1E-3),1000,1E-3);
		double right = bisect(max,1,ymax+FastMath.log(1E-3),1000,1E-3);
		
		
		double[] x = new double[100];
		for (int i=0; i<x.length; i++)
			x[i] = left+i/99.*(right-left);
		
		
		double[] fs = new double[x.length];
		double[] fs2 = new double[x.length];
		double E = 0;
		
		for (int i=0; i<fs.length; i++) 
			fs[i] = loglik(x[i])-FastMath.log(x.length-1);
		
		double[] p = fs.clone();
		
		fs[0]-=FastMath.log(2);
		fs2[0] = fs[0];
		
		for (int i=1; i<fs.length; i++) { 
			fs2[i] =  lse(fs[i-1],fs[i]-FastMath.log(2));
			fs[i] = lse(fs[i-1],fs[i]);
		}
		
		for (int i=0; i<fs.length; i++) 
			E+=x[i]*FastMath.exp(p[i]-fs2[fs2.length-1]);
		
		// normalizing
		for (int i=0; i<fs.length; i++) 
			fs2[i] -= fs2[fs.length-1];
		
		for (int i=0; i<fs.length; i++)
			fs2[i] = FastMath.exp(fs2[i]);
		
		if (betaApprox) {
			// fit beta distribution
			MultivariateFunction fun = a->{
				double s = 0;
				for (int i=0; i<x.length; i++) {
					double ss = Beta.cumulative(x[i], a[0], a[1], true, false)-fs2[i];
					s+=ss*ss;
				}
				return s;
			};
			
			// init by method of moments
			double mean = 0;
			double var = 0;
			for (int i=0; i<x.length; i++) 
				mean+=x[i]*(fs2[i]-(i==0?0:fs2[i-1]));
			for (int i=0; i<x.length; i++) {
				double d = (x[i]-mean);
				var+=d*d*(fs2[i]-(i==0?0:fs2[i-1]));
			}
			double f = (mean*(1-mean)/var-1);
			double alpha = mean*f;
			double beta = (1-mean)*f;

//				double[] theta = new SimplexOptimizer(1E-5,1E-5).optimize(
//						GoalType.MINIMIZE,
//						new NelderMeadSimplex(2),
//						new ObjectiveFunction(fun),
//						new InitialGuess(new double[]{this.alpha,this.beta}),
//						new MaxEval(10000)
//						).getPointRef();
//				
//				if (theta[0]<=0 || theta[1]<=0)
//					return null;
//			alpha = theta[0];
//			beta = theta[1];
				
				double lowerCI = Beta.quantile(lower, alpha, beta, true, false);
				double map = max;
				double upperCI = Beta.quantile(upper, alpha,beta, true, false);
				
				return new BiMixtureModelResult(lowerCI, map, upperCI, alpha, beta);
		}

		int l = Arrays.binarySearch(fs2, lower);
		if (l<0) l=Math.max(0, -l-2);
		int u = Arrays.binarySearch(fs2, upper);
		if (u<0) u=-u-1;

		
		double lowerCI = x[l];
		double map = max;
		double upperCI = x[u];
		
		return new BiMixtureModelResult(lowerCI, map, upperCI, Double.NaN, Double.NaN);
	}
	
	public final double computeMAP() {
		UnivariateFunction ufun = a->loglik(a);
		double max = new BrentOptimizer(1E-5,1E-5).optimize(
				GoalType.MAXIMIZE,
				new SearchInterval(0, 1),
				new UnivariateObjectiveFunction(ufun),
				new MaxEval(10000)
				).getPoint();
		return max;
	}
	
	private double bisect(double min, double max, double logval, int n, double eps) {
        double m;
        double fm;
        double fmin;

        fmin = loglik(min)-logval;
        
        while (true) {
        	m = UnivariateSolverUtils.midpoint(min, max);
            fm = loglik(m)-logval;

            if (fm * fmin > 0) {
                // max and m bracket the root.
                min = m;
                fmin = fm;
            } else {
                // min and m bracket the root.
                max = m;
            }

            if (FastMath.abs(max - min) <= eps) {
                m = UnivariateSolverUtils.midpoint(min, max);
                return m;
            }
        }
	}




	public final double loglik(final double x) {
		double s = 0;
		double lx = FastMath.log(x);
		double l1px = FastMath.log(1-x);
		
		
		for (BiMixtureModelData d : data) 
			s+=lse(lx+d.ll0,l1px+d.ll1)*d.count;
		
		return s;//+Beta.density(x, priorAlpha, priorBeta, true);
	}
	
	/**
	 * The elki version is much faster!
	 * @param u
	 * @param v
	 * @return
	 */
	private static final double lse(double u, double v) {
		double max = Math.max(u, v);
		final double cutoff = max - 35.350506209; // log_e(2**51)
	    double acc = 0.;
	      if(v > cutoff) {
	        acc += v < max ? FastMath.exp(v - max) : 1.;
	      }
	      if(u > cutoff) {
		        acc += u < max ? FastMath.exp(u - max) : 1.;
		      }
	    return acc > 1. ? (max + FastMath.log(acc)) : max;
	    
//		if (u>v)
//			return log1p(exp(v-u))+u;
//		else
//			return log1p(exp(u-v))+v;
	}
	
	
	
}
