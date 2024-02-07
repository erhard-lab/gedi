package gedi.slam;

import static java.lang.Math.*;

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
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.SimpleBounds;

import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.math.function.StepFunction;
import gedi.util.math.stat.NumericSample;
import gedi.util.math.stat.RandomNumbers;
import jdistlib.Beta;
import jdistlib.Binomial;

public class MCMCIntegrationProportion {

	
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
	public MCMCIntegrationProportion(double alpha, double beta, double p, double q, double lower, double upper, boolean fitbeta) {
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
		
		double[] m1 = new double[d.length];
		double[] m2 = new double[d.length];
		for (int i=0; i<m1.length; i++) {
			m1[i] = Binomial.density(d[i], n[i], q, true);
			m2[i] = Binomial.density(d[i], n[i], p, true);
		}
		double x = 0.5;
		RandomNumbers rnd = RandomNumbers.getGlobal();
		double ll = loglikprior(x,m, m1,m2);
		
		NumericSample ss = new NumericSample();
		
		for (int i=0; i<10000; i++) {
			double x1 = max(0,min(1.0,rnd.getNormal(x, 0.01)));
			double ll1 = loglikprior(x1,m, m1,m2);
			double a = Math.exp(ll1-ll);
			if (rnd.getUnif()<a) {
				x=x1;
				ll=ll1;
			}
			
			if (i>100)
				ss.add(x);
			
		}
		
		StepFunction ecdf = ss.ecdf();
		
		int l = Arrays.binarySearch(ecdf.getY(), lower);
		if (l<0) l=Math.max(0, -l-2);
		int u = Arrays.binarySearch(ecdf.getY(), upper);
		if (u<0) u=-u-1;
		
		double E = ss.evaluate(new Mean());
		
		if (fitbeta)
			return new double[] {ecdf.getX()[l],Math.min(1, Math.max(0, E)),ecdf.getX()[u],1,1};
		
		return new double[] {ecdf.getX()[l],Math.min(1, Math.max(0, E)),ecdf.getX()[u]};
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
	private double loglikprior(double x, double[] m, double[] m1, double[] m2) {
		double s = 0;
		for (int i=0; i<m1.length; i++) {
			double s1 = log(x)+m1[i];
			double s2 = log1p(-x)+m2[i];
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
