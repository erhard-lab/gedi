package gedi.slam;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.log1p;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;
import java.util.logging.Logger;

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

import gedi.util.StringUtils;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.math.stat.RandomNumbers;
import jdistlib.Beta;
import jdistlib.Binomial;

public class OptimNumericalIntegrationProportion {

	
	private double alpha;
	private double beta;
	
	private double p;
	private double q;
	private double dp;
	private double dq;
	
	private double lower;
	private double upper;
	
	private boolean fitbeta;
	
	private static final int PRE = 100;
	
	private double[][] ptab = new double[PRE][PRE];
	private double[][] qtab = new double[PRE][PRE];
	private double[][] dptab = new double[PRE][PRE];
	private double[][] dqtab = new double[PRE][PRE];
	
	
	/**
	 * 
	 * @param alpha Prior
	 * @param beta Prior
	 * @param p Mismatch prob
	 * @param q Conversion prob
	 */
	public OptimNumericalIntegrationProportion(double alpha, double beta, double p, double q, double dp, double dq, double lower, double upper, boolean fitbeta) {
		this.alpha = alpha;
		this.beta = beta;
		this.p = p;
		this.q = q;
		this.dp = dp;
		this.dq = dq;
		this.lower = lower;
		this.upper = upper;
		this.fitbeta = fitbeta;
		precompute();
	}
	
	public OptimNumericalIntegrationProportion biasedEstimator(double fp, double fq) {
		return new OptimNumericalIntegrationProportion(alpha, beta, p*fp, q*fq, dp*fp, dq*fq, lower, upper, fitbeta);
	}
	
	private void precompute() {
		// much more efficient to tabulate those all!
		for (int d=0; d<PRE; d++)
			for (int n=0; n<PRE; n++) {
				ptab[d][n] = Binomial.density(d, n, p, true); 
				qtab[d][n] = Binomial.density(d, n, q, true);
				dptab[d][n] = Binomial.density(d, n, dp, true); 
				dqtab[d][n] = Binomial.density(d, n, dq, true);
			}
	}
	
	public boolean useSingle() {
		return q>0;
	}
	public boolean useDouble() {
		return dq>0;
	}

	public SlamEstimationResult infer(double[][] x, Logger log) {
		IntArrayList d = new IntArrayList();
		IntArrayList n = new IntArrayList();
		DoubleArrayList m = new DoubleArrayList();
		
		for (int di=0; di<x.length; di++)
			for (int ni=0; ni<x.length; ni++)
				if (x[di][ni]>0) {
					n.add(ni);
					d.add(di);
					m.add(x[di][ni]);
				}
		
		return infer(d.toIntArray(),n.toIntArray(),m.toDoubleArray(),new int[0],new int[0],new double[0],log);
		
	}
	public SlamEstimationResult infer(int[] d1, int[] n, double[] m, int[] dd1, int[] dn, double[] dm, Logger log) {
		
		if (!useSingle()) d1 = new int[0];
		if (!useDouble()) dd1 = new int[0];
		
		if (d1.length==0 && dd1.length==0) {
			return new SlamEstimationResult();
		}
		int[] d = d1;
		int[] dd = dd1;
		
//		double[] m1 = new double[d.length];
//		double[] m2 = new double[d.length];
//		for (int i=0; i<m1.length; i++) {
//			m1[i] = Binomial.density(d[i], n[i], q, false);
//			m2[i] = Binomial.density(d[i], n[i], p, false);
//		}
		
//		DoubleUnaryOperator f = p->{
//			double y = 0;
//			for (int i=0; i<m1.length; i++) {
//				double r = m2[i]/(m1[i]-m2[i]);
//				y+=1/(p+r)*m[i];
//			}
//			return y;
//		};
//		DoubleUnaryOperator fprior = p->{
//			double y = f.applyAsDouble(p);
//			if (alpha!=1)
//				y+=(alpha-1)*Math.pow(p,alpha-2);
//			if (beta!=1)
//				y-=(beta-1)*Math.pow(1-p,beta-2);
//			return y;
//		};
//	
//		DoubleUnaryOperator fprime = p->{
//			double y = 0;
//			for (int i=0; i<m1.length; i++) {
//				double r = m2[i]/(m1[i]-m2[i]);
//				double pr = (p+r);
//				y+=1/(pr*pr)*m[i];
//			}
//			return -y;
//		};
//		DoubleUnaryOperator fpriorprime = p->{
//			double y = fprime.applyAsDouble(p);
//			if (alpha!=1 && alpha!=2)
//				y+=(alpha-1)*(alpha-2)*Math.pow(p,alpha-3);
//			if (beta!=1 && beta!=2)
//				y+=(beta-1)*(beta-2)*Math.pow(1-p,beta-3);
//			return y;
//		};
//
//		double max = newton(f,fprime,0.5,1000,1E-6);
//		max = newton(fprior,fpriorprime,max,1000,1E-6);
		
		// fit beta distribution
		double max = computeMAP(d,n,m,dd,dn,dm);
		
		double ymax = loglikprior(max, d, n, m,dd,dn,dm);
		
		double left = bisect(x->loglikprior(x,d,n,m,dd,dn,dm),0,max,ymax+log(1E-3),1000,1E-6);
		double right = bisect(x->loglikprior(x,d,n,m,dd,dn,dm),max,1,ymax+log(1E-3),1000,1E-6);
		
		
		double[] x = new double[100];
		for (int i=0; i<x.length; i++)
			x[i] = left+i/99.*(right-left);
		
		
		double[] fs = new double[x.length];
		double[] fs2 = new double[x.length];
		double E = 0;
		
		for (int i=0; i<fs.length; i++) 
			fs[i] = loglikprior(x[i], d, n, m,dd,dn,dm)-log(x.length-1);
		
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
			
//			try {
				double[] beta = new SimplexOptimizer(1E-5,1E-5).optimize(
						GoalType.MINIMIZE,
						new NelderMeadSimplex(2),
	//					new SimpleBounds(new double[]{1,1},new double[]{Double.POSITIVE_INFINITY,Double.POSITIVE_INFINITY}),
						new ObjectiveFunction(fun),
						new InitialGuess(new double[]{1,1}),
						new MaxEval(10000)
						).getPointRef();
				
				if (beta[0]<=0 || beta[1]<=0 && log!=null)
					log.warning("Could not fit beta distribution!");
				
				return new SlamEstimationResult(Beta.quantile(lower, beta[0], beta[1], true, false),Math.min(1, Math.max(0, E)),max,Beta.quantile(upper, beta[0], beta[1], true, false),beta[0],beta[1]);
//			} catch (TooManyEvaluationsException e) {
//				e.printStackTrace();
//				for (int i=0; i<x.length; i++) {
//					System.out.println(x[i]+"\t"+fs2[i]);
//				}
//			}
		}
		
		return new SlamEstimationResult(x[l],Math.min(1, Math.max(0, E)),max,x[u]);
	}
	
	
	public double computeMAP(int[] d, int[] n, double[] m,int[] dd, int[] dn, double[] dm) {
		UnivariateFunction ufun = a->loglikprior(a, d, n, m,dd, dn, dm);
		double max = new BrentOptimizer(1E-5,1E-5).optimize(
				GoalType.MAXIMIZE,
				new SearchInterval(0, 1),
				new UnivariateObjectiveFunction(ufun),
				new MaxEval(10000)
				).getPoint();
		return max;
	}

	private double newton(DoubleUnaryOperator f, DoubleUnaryOperator fprime, double start, int n, double eps) {
		double x0 = start;
        double x1;
        while (true) {
        	double fx = f.applyAsDouble(x0);
        	double fpx = fprime.applyAsDouble(x0);
            
            x1 = x0 - (fx / fpx);
            x1 = Math.max(0, Math.min(1,x1));
            if (FastMath.abs(x1 - x0) <= eps) {
                return x1;
            }

            x0 = x1;
            if (--n==0)
            	throw new RuntimeException("Did not converge! "+x1+" "+x0+" "+f.applyAsDouble(start));
        }
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


	private double loglikprior(final double x, int[] d, int[] n, double[] m, int[] dd, int[] dn, double[] dm) {
		double s = 0;
		double lx = log(x);
		double l1px = log1p(-x);
		for (int i=0; i<d.length; i++) {
			if (n[i]<PRE) {
				double s1 = lx+qtab[d[i]][n[i]]; //Binomial.density(d[i], n[i], q, true);
				double s2 = l1px+ptab[d[i]][n[i]]; //Binomial.density(d[i], n[i], p, true);
				s+=lse(s1,s2)*m[i];	
			}
			else {
				double s1 = lx+Binomial.density(d[i], n[i], q, true);
				double s2 = l1px+Binomial.density(d[i], n[i], p, true);
				s+=lse(s1,s2)*m[i];
			}
		}
		for (int i=0; i<dd.length; i++) {
			if (dn[i]<PRE) {
				double s1 = lx+dqtab[dd[i]][dn[i]]; //Binomial.density(dd[i], dn[i], dq, true);
				double s2 = l1px+dptab[dd[i]][dn[i]]; //Binomial.density(dd[i], dn[i], dp, true);
				s+=lse(s1,s2)*dm[i];
			}
			else {
				double s1 = lx+Binomial.density(dd[i], dn[i], dq, true);
				double s2 = l1px+Binomial.density(dd[i], dn[i], dp, true);
				s+=lse(s1,s2)*dm[i];
			}
		}
		return s+Beta.density(x, alpha, beta, true);
	}
	
	
	private static final double lse(double u, double v) {
		if (u>v)
			return log1p(exp(v-u))+u;
		else
			return log1p(exp(u-v))+v;
	}
	
	
	public static void main(String[] args) {
		
		double p = 0.001;
		double q = 0.02;
		double mm = 0.01;
		int[] n = RandomNumbers.getGlobal().getNormal(100, 40, 3).mapToInt(d->(int)Math.round(d)).toIntArray();
		double[] mix = RandomNumbers.getGlobal().getUnif(n.length).toDoubleArray();
		int[] d = new int[n.length];
		double[] m = RandomNumbers.getGlobal().getBinom(n.length,10, 0.2).toDoubleArray();
		
		for (int i=0; i<m.length; i++)
			d[i] = RandomNumbers.getGlobal().getBinom(n[i], mix[i]<mm?q:p);
		
		
		System.out.println(StringUtils.toString(new OptimNumericalIntegrationProportion(1, 1, p, q, 0,0, 0.05, 0.95, true).infer(d, n, m, new int[0],new int[0],new double[0],null)));
		System.out.println(StringUtils.toString(new NumericalIntegrationProportion(1, 1, p, q, 0.05, 0.95, true).infer(d, n, m, null)));
		System.out.println(StringUtils.toString(new MCMCIntegrationProportion(1, 1, p, q, 0.05, 0.95, true).infer(d, n, m, null)));
		System.out.println(StringUtils.toString(new VariationalBayesProportion(1, 1, p, q).infer(d, n, m, null)));
		
	}
	
}
