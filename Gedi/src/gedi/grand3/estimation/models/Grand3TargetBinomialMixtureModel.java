package gedi.grand3.estimation.models;

import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

@Deprecated
public class Grand3TargetBinomialMixtureModel implements MaximumLikelihoodModel<KNMatrixElement[][]>{

	
	private Grand3BinomialMixtureModel tbbinom = new Grand3BinomialMixtureModel();
	protected double[][] buf;
	
	public Grand3TargetBinomialMixtureModel(int numSub) {
		buf = new double[numSub][3];
	}
	
	public Grand3TargetBinomialMixtureModel set(ModelStructure[] sub) {
		for (int i=0; i<sub.length; i++)
			buf[i] = sub[i].getBinom().getParameter();
		return this;
	}
	
	public boolean set(ModelStructure[][][] models, int t, int i) {
		// either all or no model exists (i.e. the label is not there for the condition ,e.g. no4sU)
		int present = 0;
		for (int s=0; s<buf.length; s++) {
			if (models[s][t][i]!=null) {
				System.arraycopy(models[s][t][i].getBinom().getParameter(), 0, buf[s], 0, buf[s].length);
				present++;
			}
		}
		if (present>0 && present!=buf.length)
			throw new RuntimeException("Fatal error: Either all or no model should exist for subreads!");
		return present>0;
	}
	
	
	@Override
	public double logLik(KNMatrixElement[][] data, double[] param) {
		double re = 0;
		for (int i=0; i<buf.length; i++) {
			buf[i][0] = param[0];
			re+=tbbinom.logLik(data[i], buf[i]);
		}
		return re;
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
		
		return new MaximumLikelihoodParametrization("binomTargetNtr",paramNames,params,paramFormat,lowerBounds,upperBounds);
	}



	


}
