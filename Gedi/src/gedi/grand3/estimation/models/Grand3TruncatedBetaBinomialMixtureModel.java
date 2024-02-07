package gedi.grand3.estimation.models;

import static java.lang.Math.exp;
import static java.lang.Math.log;

import gedi.grand3.estimation.IncompleteBeta;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import jdistlib.rng.RandomEngine;

public class Grand3TruncatedBetaBinomialMixtureModel extends Grand3Model {
	
	
	
	public static double ldtbbinommix(int k, int n, double ntr, double perr, double pmconv, double shape) {
		return lse(lbinom(k,n,perr)+log(1-ntr),ldtbbinom(k,n,perr,pmconv,shape)+log(ntr));
	}
	
	private static double expectedTbBeta(int l, double u, double s1, double s2) {
		return exp(logdiff(IncompleteBeta.libeta(u, 1+s1, s2),IncompleteBeta.libeta(l, 1+s1, s2))-logdiff(IncompleteBeta.libeta(u, s1, s2),IncompleteBeta.libeta(l, s1, s2)));
	}

	public MaximumLikelihoodParametrization createPar() {
		return createPar(0.1,4E-4,0.05,4.5);
	}

	public MaximumLikelihoodParametrization createPar(double ntr, double perr, double pmconv, double shape) {
		return new MaximumLikelihoodParametrization(
				"tbbinom",
				new String[] {"ntr","p.err","p.mconv","shape"},
				new double[] {ntr,perr,pmconv,shape},
				new String[] {"%.4f","%.4g","%.4f","%.2f"},
				new double[] {0,0,0,-2},
				new double[] {1,0.01,0.2,5}
				);
	}


	@Override
	public double computeNewComponentExpectedProbability(MaximumLikelihoodParametrization par) {
		return expectedTbBeta(0,par.getParameter(2),exp(par.getParameter(3)),exp(-par.getParameter(3)));
	}

	@Override
	public double logLik(int k, int n, double[] param) {
		return ldtbbinommix(k, n, param[0], param[1], param[2], param[3]);
	}
	
	@Override
	public double logLikOld(int k, int n, double[] param) {
		return lbinom(k, n, param[1]);
	}
	
	@Override
	public double logLikNew(int k, int n, double[] param) {
		return ldtbbinom(k,n,param[1],param[2],param[3]);
	}
	
	@Override
	protected void sampleFromComponent1(MaximumLikelihoodParametrization par, int n, RandomEngine rnd, double[] re, int start, int end) {
		rbinom(n, par.getParameter(1),rnd, re, start, end);
	}


	@Override
	protected void sampleFromComponent2(MaximumLikelihoodParametrization par, int n, RandomEngine rnd, double[] re, int start, int end) {
		rtbbinom(n, par.getParameter(1), par.getParameter(2), par.getParameter(3),rnd, re, start, end);
	}

	

}
