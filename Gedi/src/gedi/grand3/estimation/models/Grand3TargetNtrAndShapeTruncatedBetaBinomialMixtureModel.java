package gedi.grand3.estimation.models;

import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class Grand3TargetNtrAndShapeTruncatedBetaBinomialMixtureModel extends Grand3TargetTruncatedBetaBinomialMixtureModel {

	
	
	public Grand3TargetNtrAndShapeTruncatedBetaBinomialMixtureModel(int numSub) {
		super(numSub);
	}

	@Override
	public MaximumLikelihoodParametrization createPar() {
		String[] paramNames = new String[2];
		double[] params = new double[2];
		String[] paramFormat = new String[2];
		double[] lowerBounds = new double[2];
		double[] upperBounds = new double[2];
		paramNames[0] = "ntr";
		params[0] = 0;
		paramFormat[0] = "%.4f";
		lowerBounds[0] = 0;
		upperBounds[0] = 1;
		
		paramNames[1] = "shape";
		params[1] = 0;
		paramFormat[1] = "%.2f";
		lowerBounds[1] = -2;
		upperBounds[1] = 5;
				
		for (int s=0; s<buf.length; s++) {
			params[0]+=buf[s][0];
			params[1]+=buf[s][3];
		}
		params[0]/=buf.length;
		params[1]/=buf.length;
		
		return new MaximumLikelihoodParametrization("tbbinomTargetNtrAndShape",paramNames,params,paramFormat,lowerBounds,upperBounds);
	}

	@Override
	protected void adaptBuf(double[] buf, double[] param) {
		buf[0] = param[0];
		buf[3] = param[1];
	}

	


}
