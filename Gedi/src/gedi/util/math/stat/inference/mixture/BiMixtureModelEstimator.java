package gedi.util.math.stat.inference.mixture;


import java.io.IOException;
import java.net.UnknownHostException;
import java.util.Arrays;
import java.util.function.ToDoubleFunction;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.UnivariateSolverUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.linear.LinearConstraint;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer.PopulationSize;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer.Sigma;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.random.MersenneTwister;

import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.math.optim.NelderMead;
import gedi.util.math.optim.NelderMead.NelderMeadResult;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.mutable.MutableInteger;
import gedi.util.r.RDataWriter;
import gedi.util.r.RProcess;
import jdistlib.Beta;
import jdistlib.math.MultivariableFunction;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.BobyqaConfig;
import jdistlib.math.opt.OptimizationConfig;
import jdistlib.math.opt.OptimizationResult;
import net.jafama.FastMath;

public class BiMixtureModelEstimator {

	private BiMixtureModelData[] data;
//	private double priorAlpha = 1;
//	private double priorBeta = 1;
	
	private RandomNumbers rng = RandomNumbers.getGlobal();
	private boolean fitMix = true;
	

	public BiMixtureModelEstimator(BiMixtureModelData[] data) {
		this.data = data;
	}
	
	
	public void setFitMix(boolean fitMix) {
		this.fitMix = fitMix;
	}
	
	public void setRandom(RandomNumbers rng) {
		this.rng = rng;
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
		
		double ll_0 = loglik(0);
		double ll_1 = loglik(1);
		double ymax = loglik(max);
		
		double left;
		double right;

		if (ll_0>ymax && ll_1>max) {
			
			if (ll_0>ll_1) {
				max = 0;
				ymax = ll_0;
			} else {
				max = 1;
				ymax = ll_1;
			}
			
			left = 0;
			right = 1;
			
		} else if (ll_0>ymax) {
			max = 0;
			ymax = ll_0;
			left = 0;
			right = ymax-ll_1<FastMath.log(1E3)?1:bisect(max,1,ymax+FastMath.log(1E-3),1000,1E-3);
		} else if (ll_1>ymax) {
			max = 1;
			ymax = ll_1;
			left = ymax-ll_0<FastMath.log(1E3)?0:bisect(0,max,ymax+FastMath.log(1E-3),1000,1E-3);
			right = 1;
		} else {
			left = ymax-ll_0<FastMath.log(1E3)?0:bisect(0,max,ymax+FastMath.log(1E-3),1000,1E-3);
			right = ymax-ll_1<FastMath.log(1E3)?1:bisect(max,1,ymax+FastMath.log(1E-3),1000,1E-3);
		}
		
		double[] x = new double[100];
		for (int i=0; i<x.length; i++)
			x[i] = left+i/99.*(right-left);
		
		double inte = 0;
		double[] fs = new double[x.length];
		double[] fs2 = new double[x.length];
//		double E = 0;
		
		for (int i=0; i<fs.length; i++) {
			double l = loglik(x[i]);
			fs[i] = l-FastMath.log(x.length-1);
			if (i==0)
				inte = l;
			else
				inte = lse(inte,l);
		}
		
		double[] p = fs.clone();
		
		fs[0]-=FastMath.log(2);
		fs2[0] = fs[0];
		
		for (int i=1; i<fs.length; i++) { 
			fs2[i] =  lse(fs[i-1],fs[i]-FastMath.log(2));
			fs[i] = lse(fs[i-1],fs[i]);
		}
		
//		for (int i=0; i<fs.length; i++) 
//			E+=x[i]*FastMath.exp(p[i]-fs2[fs2.length-1]);
		
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
			ToDoubleFunction<double[]> fun2 = a->{
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

			double[] theta = new NelderMead().minimize(new double[]{alpha,beta}, fun2).addLowerBound(0, 0).addLowerBound(1, 0).minimize().getX();
//			
//				double[] theta = new SimplexOptimizer(1E-5,1E-5).optimize(
//						GoalType.MINIMIZE,
//						new NelderMeadSimplex(2),
//						new ObjectiveFunction(fun),
//						new InitialGuess(new double[]{alpha,beta}),
//						new MaxEval(10000)
//						).getPointRef();
				
				alpha = theta[0];
				beta = theta[1];
				
				double lowerCI = Beta.quantile(lower, alpha, beta, true, false);
				double map = max;
				double upperCI = Beta.quantile(upper, alpha,beta, true, false);
				
				inte = inte+FastMath.log(x[1]-x[0]);
				
				double[] mix = null;
				if (fitMix) mix = fitBetaMixture();
				
				return new BiMixtureModelResult(lowerCI, map, upperCI, alpha, beta,inte,mix);
				
		}

		int l = Arrays.binarySearch(fs2, lower);
		if (l<0) l=Math.max(0, -l-2);
		int u = Arrays.binarySearch(fs2, upper);
		if (u<0) u=-u-1;

		
		double lowerCI = x[l];
		double map = max;
		double upperCI = x[u];
		
		return new BiMixtureModelResult(lowerCI, map, upperCI, Double.NaN, Double.NaN, Double.NaN,null);
	}
	
	/**
	 * 
	 * @return ntr, w, a1, b1, a2, b2, SSE, integral
	 */
	public double[] fitBetaMixture() {
		double[] ll0 = new double[data.length];
		double[] ll1 = new double[data.length];
		double[] c = new double[data.length];
		for (int i=0; i<data.length; i++) {
			ll0[i] = data[i].ll0;
			ll1[i] = data[i].ll1;
			c[i] = data[i].count;
		}
		
		
		RProcess proc = RProcess.getForCurrentThread();
		try {
			RDataWriter setter = proc.startSetting();
			
			setter.write("ll0", ll0);
			setter.write("ll1", ll1);
			setter.write("c", c);
			setter.finish();
			
			proc.library("grandR");
			double[] re = proc.callNumericFunction("grandR:::fit.ntr.betamix(ll0,ll1,c)");
			return re;
			
		} catch (Exception e) {
			throw new RuntimeException("Cannot call R!",e);
		}
		
		
	}
			
//	public double[] fitBetaMixture(double[] x, double[] Femp,
//			int nstart, double minp, double maxp, double shapeLower) {
//
//		if (x.length != Femp.length) throw new IllegalArgumentException("x and Femp must be same length.");
//		for (int i = 1; i < x.length; i++) if (x[i] <= x[i - 1]) throw new IllegalArgumentException("x must be strictly increasing.");
//
//		// Method-of-moments init
//		int n = x.length;
//		double[] dx = new double[n - 1];
//		double[] xm = new double[n - 1];
//		double[] dens = new double[n - 1];
//		double[] w = new double[n - 1];
//
//		for (int i = 0; i < n - 1; i++) {
//			dx[i] = x[i + 1] - x[i];
//			xm[i] = (x[i] + x[i + 1]) / 2;
//			dens[i] = Math.max((Femp[i + 1] - Femp[i]) / dx[i], 0);
//			w[i] = dens[i] * dx[i];
//		}
//		int[] ord = ArrayUtils.seq(0, xm.length-1, 1);
//		ArrayUtils.parallelSort(xm, ord);
//		double cumw = 0;
//		double cut = 0;
//		for (int i = 0; i < ord.length; i++) {
//			cumw += w[ord[i]];
//			if (cumw >= 0.5) {
//				cut = xm[ord[i]];
//				break;
//			}
//		}
//
//		int[] cl = new int[n - 1];
//		for (int i = 0; i < xm.length; i++) cl[i] = (xm[i] <= cut) ? 1 : 2;
//
//		double[][] ab_init = new double[2][2];
//		for (int j = 1; j <= 2; j++) {
//			double sumW = 0, sumX = 0, sumX2 = 0;
//			for (int i = 0; i < cl.length; i++) {
//				if (cl[i] == j) {
//					sumW += w[i];
//					sumX += xm[i] * w[i];
//				}
//			}
//			double m = sumX / sumW;
//			for (int i = 0; i < cl.length; i++) {
//				if (cl[i] == j) {
//					sumX2 += Math.pow(xm[i] - m, 2) * w[i];
//				}
//			}
//			double v = sumX2 / sumW;
//			double t = (m * (1 - m) / v) - 1;
//			ab_init[j - 1][0] = Math.max(m * t, shapeLower);
//			ab_init[j - 1][1] = Math.max((1 - m) * t, shapeLower);
//		}
//
//		double p_init = 0;
//		for (int i = 0; i < cl.length; i++) if (cl[i] == 1) p_init += w[i];
//
//		double m1 = ab_init[0][0] / (ab_init[0][0] + ab_init[0][1]);
//		double m2 = ab_init[1][0] / (ab_init[1][0] + ab_init[1][1]);
//		if (m1 > m2) {
//			double[] tmp = ab_init[0];
//			ab_init[0] = ab_init[1];
//			ab_init[1] = tmp;
//			p_init = 1 - p_init;
//		}
//double eps= 10E-4;
//		double[] mom_init = new double[]{p_init, ab_init[0][0], ab_init[0][1], ab_init[1][0], ab_init[1][1]};
//		double[] bestParams = null;
//		double bestSSE = Double.POSITIVE_INFINITY;
//		
//		MutableInteger in = new MutableInteger();
//
//		ToDoubleFunction<double[]> sseFunction = par->{
//			double w1 = par[0], a1 = par[1], b1 = par[2], a2 = par[3], b2 = par[4];
//			double sse = 0;
//			for (int j = 0; j < x.length; j++) {
//				double Fmix = w1 * Beta.cumulative(x[j],a1,b1,true,false) + (1 - w1) * Beta.cumulative(x[j],a2,b2,true,false);
//				sse += Math.pow(Femp[j] - Fmix, 2);
//			}
//			in.N++;
//			return sse;
//		};
//		
//		for (int i = 0; i < nstart; i++) {
//			double[] init;
//			if (i == 0) {
//				init = mom_init.clone();
//			} else if (i < nstart * 0.6) {
//				init = mom_init.clone();
//				for (int j = 0; j < init.length; j++) {
//					double jitter = Math.max(0.1, 0.2 * Math.abs(init[j]));
//					init[j] += rng.getNormal() * jitter;
//					init[j] = Math.min(Math.max(init[j], j == 0 ? minp*(1+eps) : shapeLower*(1+eps)), j == 0 ? maxp*(1-eps) : 50);
//				}
//			} else {
//				init = new double[]{
//						minp + rng.getUnif() * (maxp - minp),
//						shapeLower + rng.getUnif() * (50 - shapeLower),
//						shapeLower + rng.getUnif() * (50 - shapeLower),
//						shapeLower + rng.getUnif() * (50 - shapeLower),
//						shapeLower + rng.getUnif() * (50 - shapeLower)
//				};
//				double mean1 = init[1] / (init[1] + init[2]);
//				double mean2 = init[3] / (init[3] + init[4]);
//				if (mean1 > mean2) {
//					double tmp1 = init[1], tmp2 = init[2];
//					init[1] = init[3];
//					init[2] = init[4];
//					init[3] = tmp1;
//					init[4] = tmp2;
//				}
//			}
//
//			PointValuePair pvp = new BOBYQAOptimizer(11).optimize(
//					new MaxEval(50000),
//		            GoalType.MINIMIZE,
//		            new InitialGuess(init),
//		            new ObjectiveFunction(new MultivariateFunction() {
//						
//						@Override
//						public double value(double[] par) {
//							double w1 = par[0], a1 = par[1], b1 = par[2], a2 = par[3], b2 = par[4];
//							double sse = 0;
//							for (int j = 0; j < x.length; j++) {
//								double Fmix = w1 * Beta.cumulative(x[j],a1,b1,true,false) + (1 - w1) * Beta.cumulative(x[j],a2,b2,true,false);
//								double x = Femp[j] - Fmix;
//								sse += x*x;
//							}
//							in.N++;
//							return sse;
//						}
//					}),
//		            new SimpleBounds(new double[] {minp,shapeLower,shapeLower,shapeLower,shapeLower},
//					new double[] {maxp,9999,9999,9999,9999})
//					);
//			
////			PointValuePair pvp = new CMAESOptimizer(50000,1e-9,true,0,0,new MersenneTwister(),false,null).optimize(
////					new MaxEval(50000),
////		            GoalType.MINIMIZE,
////		            new InitialGuess(init),
////		            new ObjectiveFunction(new MultivariateFunction() {
////						
////						@Override
////						public double value(double[] par) {
////							double w1 = par[0], a1 = par[1], b1 = par[2], a2 = par[3], b2 = par[4];
////							double sse = 0;
////							for (int j = 0; j < x.length; j++) {
////								double Fmix = w1 * Beta.cumulative(x[j],a1,b1,true,false) + (1 - w1) * Beta.cumulative(x[j],a2,b2,true,false);
////								double x = Femp[j] - Fmix;
////								sse += x*x;
////							}
////							in.N++;
////							return sse;
////						}
////					}),
////		            new SimpleBounds(new double[] {minp,shapeLower,shapeLower,shapeLower,shapeLower},
////					new double[] {maxp,9999,9999,9999,9999}),
////		            new PopulationSize(10),     // must be >= 4
////		            new Sigma(new double[]{0.05, 0.2,0.2,0.2,0.2})
////					);
//			
//			if (pvp.getValue() < bestSSE) {
//				bestSSE = pvp.getValue();
//				bestParams = pvp.getPoint();
//			}
//			
////	OptimizationResult opt = Bobyqa.bobyqa(init, 
////					new double[] {minp,shapeLower,shapeLower,shapeLower,shapeLower},
////					new double[] {maxp,9999,9999,9999,9999},
////					new MultivariableFunction() {
////						@Override
////						public double eval(double... par) {
////							double w1 = par[0], a1 = par[1], b1 = par[2], a2 = par[3], b2 = par[4];
////							double sse = 0;
////							for (int j = 0; j < x.length; j++) {
////								double Fmix = w1 * Beta.cumulative(x[j],a1,b1,true,false) + (1 - w1) * Beta.cumulative(x[j],a2,b2,true,false);
////								sse += Math.pow(Femp[j] - Fmix, 2);
////							}
////							in.N++;
////							return sse;
////						}}, init.length + 6, 0.1,
////					1e-7, 50000, true);
////			if (opt.mF < bestSSE) {
////				bestSSE = opt.mF;
////				bestParams = opt.mX;
////			}
//			
////			NelderMeadResult result = new NelderMead().minimize(init, sseFunction)
////					.addLowerBound(0,minp)
////					.addLowerBound(1,shapeLower)
////					.addLowerBound(2,shapeLower)
////					.addLowerBound(3,shapeLower)
////					.addLowerBound(4,shapeLower)
////					.addUpperBound(0, maxp).minimize();
//
////			if (result.getY() < bestSSE) {
////				bestSSE = result.getY();
////				bestParams = result.getX();
////			}
//		}
//
//		if (bestParams == null) throw new RuntimeException("All starts failed");
//		System.out.println(Arrays.toString(data));
//		System.out.println(Arrays.toString(bestParams));
//		System.out.println(in);
//		
//		return bestParams;
//	}
	
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
	
	
	public static void main(String[] args) {
//		BinomialDistribution e10 = new BinomialDistribution(10, 4E-4);
//		BinomialDistribution e11 = new BinomialDistribution(11, 4E-4);
		BinomialDistribution m10 = new BinomialDistribution(10, 0.05);
		BinomialDistribution m11 = new BinomialDistribution(11, 0.05);
//		//> set.seed(42)
//		//> mm=CreateMixMatrix(n.vector=c(0,0,0,0,0,0,0,0,0,10,10),par=model.par(ntr=0.5))
//		//> fit.ntr(mm,model.par(),beta.approx = TRUE,plot=TRUE)
		BiMixtureModelData[] data = {
				new BiMixtureModelData(new BinomialDistribution(8, 0.0473).logProbability(0), new BinomialDistribution(8, 0.003369).logProbability(0),1),
				new BiMixtureModelData(new BinomialDistribution(9, 0.0473).logProbability(0), new BinomialDistribution(9, 0.003369).logProbability(0),4),
				new BiMixtureModelData(new BinomialDistribution(15, 0.0473).logProbability(0), new BinomialDistribution(15, 0.003369).logProbability(0),3),
				new BiMixtureModelData(new BinomialDistribution(16, 0.0473).logProbability(0), new BinomialDistribution(16, 0.003369).logProbability(0),1),
				new BiMixtureModelData(new BinomialDistribution(18, 0.0473).logProbability(0), new BinomialDistribution(18, 0.003369).logProbability(0),3),
				new BiMixtureModelData(new BinomialDistribution(20, 0.0473).logProbability(0), new BinomialDistribution(20, 0.003369).logProbability(0),1)
		};
		for (BiMixtureModelData d : data)
			System.out.println(d);
		
		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
		
//		BinomialDistribution e = new BinomialDistribution(21, 0.0001);
//		BinomialDistribution m = new BinomialDistribution(21, 0.04);
//		BiMixtureModelData[] data = {
//				new BiMixtureModelData(m.logProbability(0), e.logProbability(0), 1)
//		};
		

		
		
		
		BiMixtureModelData[] data2 = {
				new BiMixtureModelData(-0.43636748110048107,-0.030375366392205783, 1),
				new BiMixtureModelData(-0.304628430698524,-0.0037972450581003035, 1),
				new BiMixtureModelData(-0.34814677794117027,-0.004339708637828918, 1),
				new BiMixtureModelData(-1.0009219865808647,-0.012476662333758139, 1),
				new BiMixtureModelData(-0.6527752086396943,-0.008136953695929221, 1),
				new BiMixtureModelData(-0.6962935558823405,-0.008679417275657835, 1)
		};
		System.out.println(new BiMixtureModelEstimator(data2).estimate(0.95, true));
//
//		long s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
//		s = System.nanoTime();
//		System.out.println(new BiMixtureModelEstimator(data).estimate(0.95, true));
//		System.out.println(StringUtils.getHumanReadableTimespanNano(System.nanoTime()-s));
	}
}
