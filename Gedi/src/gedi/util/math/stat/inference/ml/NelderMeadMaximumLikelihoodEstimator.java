package gedi.util.math.stat.inference.ml;

import java.io.IOException;
import java.util.Vector;
import java.util.concurrent.locks.ReentrantLock;

import gedi.util.functions.EI;
import gedi.util.io.text.LineWriter;
import gedi.util.math.optim.NelderMead;
import gedi.util.math.optim.NelderMead.NelderMeadResult;
import gedi.util.math.optim.NelderMead.NelderMeadSetup;
import gedi.util.mutable.MutablePair;
import jdistlib.ChiSquare;


/**
 * This will use the nelder mead algorithm with a 1E-4 absolute convergence criterion (i.e. if nm cannot improve the current value by a factor of 1E-4m it will stop)
 * @author flo
 *
 */
public class NelderMeadMaximumLikelihoodEstimator implements MaximumLikelihoodEstimator {
	

	
	private double confidenceLevel = 0.975;
	private int nthreads = Math.max(Runtime.getRuntime().availableProcessors(),8);
	private double bounds_eps = 1E-4;
	private int profilePoints = 100;
	
	
	public NelderMeadMaximumLikelihoodEstimator() {
	}
	public NelderMeadMaximumLikelihoodEstimator(double confidenceLevel, int nthreads, double bounds_eps, int profilePoints) {
		this.confidenceLevel = 1-0.5*(1-confidenceLevel);
		this.nthreads = nthreads;
		this.bounds_eps = bounds_eps;
		this.profilePoints = profilePoints;
	}
	
	
	public void setNthreads(int nthreads) {
		this.nthreads = nthreads;
	}

	protected NelderMeadSetup addParameterBounds(MaximumLikelihoodParametrization par, NelderMeadSetup setup, boolean[] fixed) {
		for (int i=0; i<par.params.length; i++) {
			if (fixed==null || !fixed[i]) {
				if (Double.isFinite(par.lowerBounds[i])) setup.addLowerBound(i, par.lowerBounds[i]);
				if (Double.isFinite(par.upperBounds[i])) setup.addUpperBound(i, par.upperBounds[i]);
			}
		}
		return setup;
	}
	
	@Override
	public <D> MaximumLikelihoodParametrization fit(MaximumLikelihoodModel<D> model, MaximumLikelihoodParametrization par, D data, String... fixedNames) {
		for (int p=0; p<par.params.length; p++) 
			par.params[p] = Math.max(par.lowerBounds[p]*(1+bounds_eps), Math.min(par.upperBounds[p]*(1-bounds_eps),par.params[p]));
		
		boolean[] fixed = new boolean[par.params.length];
		for (String n : fixedNames)
			fixed[par.getParameterIndex(n)]=true;
		
		
		NelderMead nm = new NelderMead();
//		nm.setReltol(0);
		nm.setAbstol(1E-4);
		nm.setFnscale(-1);
		nm.setMaxit(1000);
		
		double[] param = par.params.clone();
		double start = model.logLik(data,param);
		
		NelderMeadResult res = addParameterBounds(par,nm.minimize(
				param, 
				a->model.logLik(data,a)-start), fixed)
				.setFixedParameter(fixed)
				.minimize();
		
        double[] mc = res.getX();
        par.logLik = res.getY()+start;
        System.arraycopy(mc, 0, par.params, 0, param.length);
        par.convergence = res.getConvergence();
        	
        return par;
	}
	
	/**
	 * Write prefix in front of each line (e.g. "condition\t")
	 * @param out
	 * @param prefix
	 * @param points
	 * @param nthreads
	 * @return
	 * @throws IOException
	 */
	@Override
	public <D> MaximumLikelihoodParametrization writeProfileLikelihoods(MaximumLikelihoodModel<D> model, MaximumLikelihoodParametrization par, D  data, LineWriter out, String prefix, String... fixedNames) throws IOException {
		par.resetConvergence();
		
		if (!par.isInBounds()) 
			return par;
		
		boolean[] fixed = new boolean[par.params.length];
		for (String n : fixedNames)
			fixed[par.getParameterIndex(n)]=true;

		if (par.lowerBounds==null) throw new RuntimeException("Call computeProfileLikelihoodPointwiseConfidenceIntervals first!");
		
		NelderMead nm = new NelderMead();
//		nm.setReltol(0);
		nm.setAbstol(1E-4);
		nm.setFnscale(-1);
		nm.setMaxit(1000);
		double opt = model.logLik(data, par);//par.logLik;
		Vector<String> errs = new Vector<>();
		
		EI.seq(0, par.params.length).filter(p->!fixed[p]).unfold(p->{
			double l = Math.max(par.lowerBounds[p]*(1+1E-4), par.params[p]-(par.params[p]-par.lowerConfidence[p])*1.1);
			double u = Math.min(par.upperBounds[p]*(1-1E-4), par.params[p]+(par.upperConfidence[p]-par.params[p])*1.1);
			double step = (u-l)/profilePoints;
			return EI.seq(l, u+step/2,step).map(v->new MutablePair<>(p, v));
		}).parallelized(nthreads, 1, ei->ei.map(pair->{
			int p = pair.Item1;
			double v = pair.Item2;
			
			double[] x = par.params.clone();
			boolean[] hfixed = new boolean[par.paramNames.length];
			hfixed[p] = true;
			for (int i=0; i<fixed.length; i++)
				hfixed[i] |= fixed[i];
			x[p]=v;
			NelderMeadResult res = addParameterBounds(par,nm.minimize(
					x, a->model.logLik(data,a)-opt), hfixed)
					.setFixedParameter(hfixed)
					.minimize();
			if (res.getConvergence()!=null && res.getConvergence().length()>0) {
				par.addConvergence("Could not compute profile likelihood: "+res.getConvergence());
				errs.add(par.getConvergence());
				return null;
			}
			return prefix+par.paramNames[p]+"\t"+EI.seq(0,x.length).map(i->fixed[i]?"NA":res.getX()[i]).concat("\t")+"\t"+res.getY();
		})).removeNulls().print(out);
		
		if(errs.size()>0) {
			par.resetConvergence();
			par.addConvergence(errs.get(0));
		}
		
		return par;
	}

	
	
	@Override
	public <D> MaximumLikelihoodParametrization computeProfileLikelihoodPointwiseConfidenceIntervals(MaximumLikelihoodModel<D> model, MaximumLikelihoodParametrization par, D  data, String... fixedNames) {
		par.resetConvergence();
		
		if (!par.hasConfidences()) {
			par.lowerConfidence = new double[par.params.length];
			par.upperConfidence = new double[par.params.length];
		}
		
		if (par.isInBounds()) {
			boolean[] fixed = new boolean[par.params.length];
			for (String n : fixedNames)
				fixed[par.getParameterIndex(n)]=true;
	
			double llDrop = ChiSquare.quantile(confidenceLevel, 1, true, false);
			
			
			double opt = model.logLik(data,par.params);
			double[] lower = par.lowerBounds;
			double[] upper = par.upperBounds;
			ReentrantLock lock = new ReentrantLock();
			
			EI.seq(-par.params.length, par.params.length).parallelized(Math.min(nthreads, par.params.length*2), 1, ei->ei.map(index->{
				boolean islow = index<0;
				int p=islow?-index-1:index;
				if (!fixed[p]) {
					double val = islow?lineSearch(model,par,data,par.params.clone(),p,lower,fixed,-1,llDrop, opt):lineSearch(model,par,data,par.params.clone(),p,upper,fixed,1,llDrop, opt);
					lock.lock();
					if (islow) 
						par.setLowerConfidence(p, val);
					else
						par.setUpperConfidence(p, val);
					lock.unlock();
				}
				return 1;
			})).drain();
		} 
		else {
			par.setConfidenceToBounds();
		}
			
		return par;
	}
	
	
	private <D> double lineSearch(MaximumLikelihoodModel<D> model, MaximumLikelihoodParametrization par, D  data, double[] x, int p, double[] bounds, boolean[] fixed, int f, double llDrop, double opt) {
		
		NelderMead nm = new NelderMead();
//		nm.setReltol(0);
		nm.setAbstol(1E-4);
		nm.setFnscale(-1);
		nm.setMaxit(1000);
		
		double[] ox = x.clone();
		double d = 0.001*x[p];
		boolean[] hfixed = new boolean[par.paramNames.length];
		hfixed[p] = true;
		for (int i=0; i<fixed.length; i++)
			hfixed[i] |= fixed[i];
		
		for (;;) {
			x[p]=Math.min(f*bounds[p], f*(ox[p]+f*d))*f;
			boolean bisect = false;
			if (x[p]<=par.lowerBounds[p]) {
				x[p] = par.lowerBounds[p];
				bisect = true;
			}
			else if (x[p]>=par.upperBounds[p]) {
				x[p] = par.upperBounds[p];
				bisect = true;
			}
			else {
				NelderMeadResult res = addParameterBounds(par,nm.minimize(
						x, a->model.logLik(data,a)-opt), hfixed)
						.setFixedParameter(hfixed)
						.minimize();
				if (res.getConvergence()!=null && res.getConvergence().length()>0) {
					par.addConvergence("Could not compute profile likelihood: "+res.getConvergence());
					return Double.NaN;
				}
				bisect = -res.getY()>llDrop;
			}
			
			if (bisect) {
				// now we can do bisection!
				double u = ox[p];
				double l = x[p];
				for (;;) {
					x[p] = (u+l)/2;
					NelderMeadResult res = addParameterBounds(par,nm.minimize(
							x, a->model.logLik(data,a)-opt), hfixed)
							.addFixedParameter(p)
							.minimize();
					if (res.getConvergence()!=null && res.getConvergence().length()>0) {
						par.addConvergence("Could not compute profile likelihood: "+res.getConvergence());
						return Double.NaN;
					}
					
					if (Math.abs(-res.getY()-llDrop)/llDrop<1E-2)
						return x[p];
					if (-res.getY()>llDrop) l=x[p];
					else u=x[p];
					if (Math.abs(u-l)/ox[p]<1E-2)
						return (u+l)/2;
				}
			}
			if (x[p]==bounds[p]) {
				return x[p];
			}
			
			d*=2;
		}
	}
	
	

}
