package gedi.grand3.estimation.models;

import static java.lang.Math.log;

import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import jdistlib.rng.RandomEngine;

public class Grand3BinomialMixtureModel extends Grand3Model {
	
	
	
	public static double ldbinommix(int k, int n, double ntr, double perr, double pconv) {
		return lse(lbinom(k,n,perr)+log(1-ntr),lbinom(k,n,pconv)+log(ntr));
	}
	
	public MaximumLikelihoodParametrization createPar() {
		return createPar(0.1,4E-4,0.05);
	}

	public MaximumLikelihoodParametrization createPar(double ntr, double perr, double pconv) {
		return new MaximumLikelihoodParametrization(
				"binom",
				new String[] {"ntr","p.err","p.conv"},
				new double[] {ntr,perr,pconv},
				new String[] {"%.4f","%.4g","%.4f"},
				new double[] {0,0,0},
				new double[] {1,0.01,0.2}
			);
	}
	
	@Override
	public double computeNewComponentExpectedProbability(MaximumLikelihoodParametrization par) {
		return par.getParameter(2);
	}

	@Override
	public double logLik(int k, int n, double[] param) {
		return ldbinommix(k, n, param[0], param[1], param[2]);
	}
	
	@Override
	public double logLikOld(int k, int n, double[] param) {
		return lbinom(k, n, param[1]);
	}
	
	@Override
	public double logLikNew(int k, int n, double[] param) {
		return lbinom(k, n, param[2]);
	}
	
	
	@Override
	protected void sampleFromComponent1(MaximumLikelihoodParametrization par,int n, RandomEngine rnd, double[] re, int start, int end) {
		rbinom(n, par.getParameter(1),rnd, re, start, end);
	}


	@Override
	protected void sampleFromComponent2(MaximumLikelihoodParametrization par,int n, RandomEngine rnd, double[] re, int start, int end) {
		rbinom(n, par.getParameter(2),rnd, re, start, end);
	}

	
	

}
