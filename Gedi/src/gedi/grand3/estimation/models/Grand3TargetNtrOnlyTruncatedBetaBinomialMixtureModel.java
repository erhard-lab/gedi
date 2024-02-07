package gedi.grand3.estimation.models;

import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class Grand3TargetNtrOnlyTruncatedBetaBinomialMixtureModel extends Grand3TargetTruncatedBetaBinomialMixtureModel {

	
	public Grand3TargetNtrOnlyTruncatedBetaBinomialMixtureModel(int numSub) {
		super(numSub);
	}
	

	@Override
	public MaximumLikelihoodParametrization createPar() {
		String[] paramNames = new String[1];
		double[] params = new double[1];
		String[] paramFormat = new String[1];
		double[] lowerBounds = new double[1];
		double[] upperBounds = new double[1];
		paramNames[0] = "ntr";
		params[0] = 0;
		paramFormat[0] = "%.4f";
		lowerBounds[0] = 0;
		upperBounds[0] = 1;
		
		for (int s=0; s<buf.length; s++) {
			params[0]+=buf[s][0];
		}
		params[0]/=buf.length;
		
		return new MaximumLikelihoodParametrization("tbbinomTargetNtr",paramNames,params,paramFormat,lowerBounds,upperBounds);
	}


	@Override
	protected void adaptBuf(double[] buf, double[] param) {
		buf[0] = param[0];
	}

	


}
