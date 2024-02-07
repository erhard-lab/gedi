package gedi.util.math.optim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.ToDoubleBiFunction;
import java.util.function.ToDoubleFunction;

import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.collections.intcollections.IntArrayList;

/**
 * Nelder-Mead
 * 
 * @author Florian Erhard
 *
 */
public class NelderMead {

	private static final double BIG = 1.0e+35;

	private int maxit = 500;
	private double abstol = Math.sqrt(2.22E-16);//Double.NEGATIVE_INFINITY;
	private double reltol = Math.sqrt(2.22E-16);
	private double alpha = 1;
	private double beta = 0.5;
	private double gamma = 2;
	private double fnscale = 1;
	
	public int getMaxit() {
		return maxit;
	}

	public void setMaxit(int maxit) {
		this.maxit = maxit;
	}

	public double getAbstol() {
		return abstol;
	}

	public void setAbstol(double abstol) {
		this.abstol = abstol;
	}

	public double getReltol() {
		return reltol;
	}

	public void setReltol(double reltol) {
		this.reltol = reltol;
	}

	public double getAlpha() {
		return alpha;
	}

	public void setAlpha(double alpha) {
		this.alpha = alpha;
	}

	public double getBeta() {
		return beta;
	}

	public void setBeta(double beta) {
		this.beta = beta;
	}

	public double getGamma() {
		return gamma;
	}

	public void setGamma(double gamma) {
		this.gamma = gamma;
	}

	public double getFnscale() {
		return fnscale;
	}

	public void setFnscale(double fnscale) {
		this.fnscale = fnscale;
	}

	public double getMu() {
		return mu;
	}

	public void setMu(double mu) {
		this.mu = mu;
	}

	public int getOuterIterations() {
		return outerIterations;
	}

	public void setOuterIterations(int outerIterations) {
		this.outerIterations = outerIterations;
	}

	public double getOuterEps() {
		return outerEps;
	}

	public void setOuterEps(double outerEps) {
		this.outerEps = outerEps;
	}

	
	/**
	 * Fixed contains the indices (into x) of parameter that are fixed!
	 * If only one parameter to optimize: switch to Brent's method!
	 * 
	 * @param x
	 * @param fun
	 * @param fixed
	 * @return
	 */
	private NelderMeadResult minimize(double[] x, ToDoubleFunction<double[]> fun, int[] fixed) {

		if (maxit <= 0) {
			return new NelderMeadResult(x, fun.applyAsDouble(x), "");
		}
		
		int[] unfixed = null;
		double[] origx = x;
		if (fixed!=null && fixed.length>0) {
			double[] fixeds = new double[fixed.length];
			unfixed = new int[x.length-fixed.length];
			double[] nx = new double[x.length-fixed.length];
			double[] parbuffer = x.clone();
			int ui = 0;
			int fi = 0;
			for (int i=0; i<x.length; i++) {
				if (fi<fixed.length && i==fixed[fi]) 
					fixeds[fi++] = x[i];
				else {
					unfixed[ui] = i;
					nx[ui++] = x[i];
				}
			}
			origx=x;
			x=nx;
			
			ToDoubleFunction<double[]> oldfun = fun;
			int[] uunfixed = unfixed;
			fun = px->{
				for (int i=0; i<uunfixed.length; i++)
					parbuffer[uunfixed[i]]=px[i];
				return oldfun.applyAsDouble(parbuffer);
			};
		}
		
		
		
		String failedMessage = "";
		double y = 0;
		
		if (x.length==1) {
			try {
				y = doBrent(x,fun);
			} catch (TooManyEvaluationsException e) {
				failedMessage = e.getMessage();
			}
		}
		else {
			int funcount = 0, H, i, j, L = 0;
			double[][] P;
			double VH, VL, VR;
	
			
	
			P = new double[x.length + 1][x.length + 2];
			y = fun.applyAsDouble(x) / fnscale;
			if (Double.isInfinite(y))
				throw new RuntimeException("function cannot be evaluated at initial parameters");
	
			
			funcount = 1;
			double tol = reltol * (Math.abs(y) + reltol);
			P[x.length][0] = y;
			for (i = 0; i < x.length; i++)
				P[i][0] = x[i];
	
			L = 1;
	
			double size = 0;
			double step = 0;
			for (i = 0; i < x.length; i++) {
				if (0.1 * Math.abs(x[i]) > step)
					step = 0.1 * Math.abs(x[i]);
			}
			if (step == 0.0)
				step = 0.1;
			for (j = 2; j <= x.length + 1; j++) {
				for (i = 0; i < x.length; i++)
					P[i][j - 1] = x[i];
	
				double nstep = step;
				while (P[j - 2][j - 1] == x[j - 2]) {
					P[j - 2][j - 1] = x[j - 2] + nstep;
					nstep *= 10;
				}
				size += nstep;
			}
	
			boolean computeVertex = true;
			do {
				if (computeVertex) {
					for (j = 0; j < x.length + 1; j++) {
						if (j + 1 != L) {
							for (i = 0; i < x.length; i++)
								x[i] = P[i][j];
							y = fun.applyAsDouble(x) / fnscale;
							if (!Double.isFinite(y))
								y = BIG;
							funcount++;
							P[x.length][j] = y;
						}
					}
					computeVertex = false;
				}
	
				VL = P[x.length][L - 1];
				VH = VL;
				H = L;
	
				for (j = 1; j <= x.length + 1; j++) {
					if (j != L) {
						y = P[x.length][j - 1];
						if (y < VL) {
							L = j;
							VL = y;
						}
						if (y > VH) {
							H = j;
							VH = y;
						}
					}
				}
//				if (VH <= VL + tol || VL <= abstol) 
//					break;
				if (VH <= VL + tol || VH <= VL + abstol) 
					break;
	
				for (i = 0; i < x.length; i++) {
					double s = -P[i][H - 1];
					for (j = 0; j < x.length + 1; j++)
						s += P[i][j];
					P[i][x.length + 1] = s / x.length;
				}
				for (i = 0; i < x.length; i++)
					x[i] = (1.0 + alpha) * P[i][x.length + 1] - alpha * P[i][H - 1];
				y = fun.applyAsDouble(x) / fnscale;
				if (!Double.isFinite(y))
					y = BIG;
				funcount++;
				VR = y;
				if (VR < VL) {
					P[x.length][x.length + 1] = y;
					for (i = 0; i < x.length; i++) {
						y = gamma * x[i] + (1 - gamma) * P[i][x.length + 1];
						P[i][x.length + 1] = x[i];
						x[i] = y;
					}
					y = fun.applyAsDouble(x) / fnscale;
					if (!Double.isFinite(y))
						y = BIG;
					funcount++;
					if (y < VR) {
						for (i = 0; i < x.length; i++)
							P[i][H - 1] = x[i];
						P[x.length][H - 1] = y;
					} else {
						for (i = 0; i < x.length; i++)
							P[i][H - 1] = P[i][x.length + 1];
						P[x.length][H - 1] = VR;
					}
				} else {
					if (VR < VH) {
						for (i = 0; i < x.length; i++)
							P[i][H - 1] = x[i];
						P[x.length][H - 1] = VR;
					}
	
					for (i = 0; i < x.length; i++)
						x[i] = (1 - beta) * P[i][H - 1] + beta * P[i][x.length + 1];
					y = fun.applyAsDouble(x) / fnscale;
					if (!Double.isFinite(y))
						y = BIG;
					funcount++;
	
					if (y < P[x.length][H - 1]) {
						for (i = 0; i < x.length; i++)
							P[i][H - 1] = x[i];
						P[x.length][H - 1] = y;
					} else {
						if (VR >= VH) {
							computeVertex = true;
							double nsize = 0;
							for (j = 0; j < x.length + 1; j++) {
								if (j + 1 != L) {
									for (i = 0; i < x.length; i++) {
										P[i][j] = beta * (P[i][j] - P[i][L - 1]) + P[i][L - 1];
										size += Math.abs(P[i][j] - P[i][L - 1]);
									}
								}
							}
							if (nsize < size) {
								size = nsize;
							} else {
								failedMessage = "Degenerate simplex";
								break;
							}
						}
					}
				}
	
			} while (funcount <= maxit);
	
			for (i = 0; i < x.length; i++)
				x[i] = P[i][L - 1];
			if (funcount > maxit)
				failedMessage = "Did not converge within " + maxit + " evals";
			y = P[x.length][L - 1] * fnscale;
		}
		
		// restore fixed parameters
		if (unfixed!=null)
			for (int i=0; i<unfixed.length; i++)
				origx[unfixed[i]]=x[i];
		
		return new NelderMeadResult(origx, y, failedMessage);
	}

	private double doBrent(double[] x, ToDoubleFunction<double[]> fun) {
		  BrentOptimizer optimizer = new BrentOptimizer(abstol, reltol);
		  double[] b = new double[1];
		  double d = Math.abs(x[0]);
		  
		  for (;;) {
			  UnivariatePointValuePair val = optimizer.optimize(new MaxEval(maxit),new UnivariateObjectiveFunction(xx->{
				  b[0]=xx; return fun.applyAsDouble(b) / fnscale;
			  }), GoalType.MINIMIZE, new SearchInterval(x[0]-d,x[0]+d));
			  if (Math.abs(x[0]-d-val.getValue())>1E-5 && Math.abs(x[0]+d-val.getValue())>1E-5) {
				  x[0] = val.getPoint();
				  return val.getValue();
			  }
			  d*=2;
		  }
	}

	private double[] sub(double[] a, double[] b, double[] re) {
		for (int i = 0; i < a.length; i++)
			re[i] = a[i] - b[i];
		return re;
	}

	private double[] mult(double[][] m, double[] a, double[] re) {
		Arrays.fill(re, 0);
		for (int i = 0; i < re.length; i++)
			for (int j = 0; j < a.length; j++)
				re[i] += a[j] * m[i][j];
		return re;
	}

	private double mu = 0.0001;
	private int outerIterations = 100;
	private double outerEps = 0.00001;
	

	private NelderMeadResult minimize(double[] x, ToDoubleFunction<double[]> fun, double[][] ui, double[] ci, int[] fixed) {

		double mu = fnscale < 0 ? -this.mu : this.mu;

		double[] uix = new double[ci.length];
		double[] gi = new double[ci.length];
		double[] giold = new double[ci.length];
		double[] tmp = new double[ci.length];

		ToDoubleBiFunction<double[], double[]> R = (rx, xold) -> {
			mult(ui, rx, uix);
			sub(uix, ci, gi);
			for (double gie : gi)
				if (gie < 0)
					return Double.NaN;
			mult(ui, xold, tmp);
			sub(tmp, ci, giold);
			double bar = 0;
			for (int i = 0; i < giold.length; i++)
				bar += giold[i] * Math.log(gi[i]) - uix[i];
			if (!Double.isFinite(bar))
				bar = Double.NEGATIVE_INFINITY;
			return fun.applyAsDouble(rx) - mu * bar;
		};

		mult(ui, x, tmp);
		sub(tmp, ci, tmp);
		for (int i=0; i<tmp.length; i++)
			if (tmp[i] <= 0 || Double.isNaN(tmp[i]))
				throw new RuntimeException("Initial value is not in the interior of the feasible region: "+i);
		double y = fun.applyAsDouble(x);
		double r = R.applyAsDouble(x, x);
		double[] thetaold = new double[x.length];
		ToDoubleFunction<double[]> bfun = rx -> R.applyAsDouble(rx, thetaold);
		double smu = Math.signum(mu);

		NelderMeadResult a = null;
		int i; double yold=0,rold=0;
		for (i = 0; i < outerIterations; i++) {
			yold = y;
			rold = r;
			System.arraycopy(x, 0, thetaold, 0, x.length);
			a = minimize(x, bfun, fixed);
			r = a.value;
			if (Double.isFinite(r) && Double.isFinite(rold) && Math.abs(r - rold) < (1e-3 + Math.abs(r)) * outerEps)
				break;
			x = a.par;
			y = fun.applyAsDouble(x);
			if (smu * y > smu * yold)
				break;
		}
		if (i == outerIterations)
			a.convergence = "Barrier algorithm ran out of iterations and did not converge";
		if (mu > 0 && y > yold)
			a.convergence = "Objective function increased at outer iteration " + i;
		if (mu < 0 && y < yold)
			a.convergence = "Objective function decreased at outer iteration " + i;
		a.value = fun.applyAsDouble(a.par);
		return a;
	}
	
	public static class NelderMeadResult {
		private double[] par;
		private double value;
		private String convergence;

		public NelderMeadResult(double[] par, double value, String convergence) {
			this.par = par;
			this.convergence = convergence;
			this.value = value;
		}
		
		public double[] getX() {
			return par;
		}
		
		public double getY() {
			return value;
		}
		public String getConvergence() {
			return convergence;
		}

		@Override
		public String toString() {
			return "f(" + Arrays.toString(par) + ")=" + value + (convergence!=null && convergence.length()>0?" "+convergence:"");
		}
	}

	
	public NelderMeadSetup minimize(double[] x, ToDoubleFunction<double[]> fun) {
		return new NelderMeadSetup(x, fun);
	}
	
	public class NelderMeadSetup {
		private IntArrayList fixed = new IntArrayList();
		private DoubleArrayList ci = new DoubleArrayList();
		private ArrayList<double[]> ui = new ArrayList<>();
		
		private double[] x;
		private ToDoubleFunction<double[]> fun;
		
		public NelderMeadSetup(double[] x, ToDoubleFunction<double[]> fun) {
			this.x = x;
			this.fun = fun;
		}

		public NelderMeadSetup addFixedParameter(int index) {
			this.fixed.add(index);
			return this;
		}
		
		public NelderMeadSetup setFixedParameter(boolean[] fixed) {
			for (int i=0; i<fixed.length; i++)
				if (fixed[i])
					addFixedParameter(i);
			return this;
		}
		
		public NelderMeadSetup addLowerBound(int index, double lower) {
			ci.add(lower);
			double[] uir = new double[x.length];
			uir[index]=1;
			ui.add(uir);
			return this;
		}
		
		public NelderMeadSetup addUpperBound(int index, double upper) {
			ci.add(-upper);
			double[] uir = new double[x.length];
			uir[index]=-1;
			ui.add(uir);
			return this;
		}
		
		public NelderMeadSetup addBounds(int index, double lower, double upper) {
			addLowerBound(index, lower);
			addUpperBound(index, upper);
			return this;
		}
		
		public NelderMeadSetup addInequality(double[] coeff, boolean isGreater, double y) {
			if (!isGreater) {
				for (int i=0; i<coeff.length; i++) coeff[i]*=-1;
				y*=-1;
			}
			ci.add(y);
			ui.add(coeff);
			return this;
		}
		
		public NelderMeadResult minimize() {
			if (ci.size()==0)
				return NelderMead.this.minimize(x, fun, fixed.toIntArray());
			return NelderMead.this.minimize(x, fun, ui.toArray(new double[0][]), ci.toDoubleArray(), fixed.toIntArray());
		}

		
	}

}
