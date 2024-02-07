package gedi.slam;

import static java.lang.Math.*;

import java.math.BigDecimal;
import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.UnivariateSolverUtils;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.util.FastMath;

import gedi.util.datastructure.collections.intcollections.IntArrayList;
import jdistlib.Beta;
import jdistlib.Binomial;
import jdistlib.math.MathFunctions;

public class LikelihoodCombinationComputer {
	
	private double[] lookup;
	
	
	public LikelihoodCombinationComputer() {
		
		lookup = new double[2000_000];
		for (int i=0; i<lookup.length; i++) {
			lookup[i] = log1p(exp(-i/(double)100_000));
		}
		
	}

	public double[] compute(double[] a, double[] b) {

		double[] re = new double[a.length+1];
		re[0] = a[0];
		re[1] = b[1];
		for (int k=1; k<a.length; k++) {
			re[k+1]=re[k]+b[k];
			for (int i=k-1; i>=1; i--) {
				re[i] = lse(re[i]+a[k],re[i-1]+b[k]);
			}
			re[0]=re[0]+a[k];
		}
		for (int i=0; i<a.length; i++)
			re[i]=re[i]-MathFunctions.lchoose(a.length, i);
		return re;
	}
	
	public double fit(double[] a, double[] b) {
		double max = ml(a, b);
		double ymax = loglikprior(max, a,b);
		
		double left = bisect(x->loglikprior(x,a,b),0,max,ymax+log(1E-3),1000,1E-6);
		double right = bisect(x->loglikprior(x,a,b),max,1,ymax+log(1E-3),1000,1E-6);
		
		
		double[] x = new double[100];
		for (int i=0; i<x.length; i++)
			x[i] = left+i/99.*(right-left);
		
		
		double[] fs = new double[x.length];
		double[] fs2 = new double[x.length];
		double E = 0;
		
		for (int i=0; i<fs.length; i++) 
			fs[i] = loglikprior(x[i], a,b)-log(x.length-1);
		
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
		
		int l = Arrays.binarySearch(fs2, 0.025);
		if (l<0) l=Math.max(0, -l-2);
		int u = Arrays.binarySearch(fs2, 0.975);
		if (u<0) u=-u-1;
		
		System.out.println(x[l]+"-"+x[u]);
	
		// fit beta distribution
		MultivariateFunction fun = ax->{
			double s = 0;
			for (int i=0; i<x.length; i++) {
				double ss = Beta.cumulative(x[i], ax[0], ax[1], true, false)-fs2[i];
				s+=ss*ss;
			}
			return s;
		};
		
//			try {
			double[] beta = new SimplexOptimizer(1E-5,1E-5).optimize(
					GoalType.MINIMIZE,
					new NelderMeadSimplex(2),
//					new SimpleBounds(new double[]{1,1},new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY}),
					new ObjectiveFunction(fun),
					new InitialGuess(new double[]{1,1}),
					new MaxEval(10000)
					).getPointRef();

				System.out.println(Arrays.toString(beta));
			
				return max;
				
//			} catch (TooManyEvaluationsException e) {
//				e.printStackTrace();
//				for (int i=0; i<x.length; i++) {
//					System.out.println(x[i]+"\t"+fs2[i]);
//				}
//			}
	}
	
	public double ml(double[] a, double[] b) {
		UnivariateFunction ufun = x->loglikprior(x, a,b);
		double max = new BrentOptimizer(1E-5,1E-5).optimize(
				GoalType.MAXIMIZE,
				new SearchInterval(0, 1),
				new UnivariateObjectiveFunction(ufun),
				new MaxEval(10000)
				).getPoint();
		return max;
	}

	private double loglikprior(final double x, double[] a, double[] b) {
		double s = 0;
		double lx = log(x);
		double l1px = log1p(-x);
		for (int i=0; i<a.length; i++) {
			double s1 = lx+a[i];
			double s2 = l1px+b[i];
			s+=lse(s1,s2);
		}
		return s;
	}
	

	private final double lse(double u, double v) {
		if (u>v)
			return calc(v-u)+u;
		else
			return calc(u-v)+v;
	}
	
	public double calc(double val) {
		if (val<=0 && val>-20) 
			return lookup[(int) (-val*100_000)];
		return log1p(exp(val));
//	    final long tmp = (long) (1512775 * val + 1072632447);
//	    return Double.longBitsToDouble(tmp << 32);
	}

	
	
	
	private double bisect(DoubleUnaryOperator f, double min, double max, double logval, int n, double eps) {
        double m;
        double fm;
        double fmin;

        while (true) {
        	m = UnivariateSolverUtils.midpoint(min, max);
            fmin = f.applyAsDouble(min)-logval;
            fm = f.applyAsDouble(m)-logval;

            if (fm * fmin > 0) {
                // max and m bracket the root.
                min = m;
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

}
