package gedi.grand3.estimation.models;

import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;

public abstract class Grand3TargetTruncatedBetaBinomialMixtureModel implements MaximumLikelihoodModel<KNMatrixElement[][]>{

	
	private Grand3TruncatedBetaBinomialMixtureModel tbbinom = new Grand3TruncatedBetaBinomialMixtureModel();
	protected double[][] buf;
	
	public Grand3TargetTruncatedBetaBinomialMixtureModel(int numSub) {
		buf = new double[numSub][4];
	}
	
	public Grand3TargetTruncatedBetaBinomialMixtureModel set(ModelStructure[] sub) {
		for (int i=0; i<sub.length; i++)
			buf[i] = sub[i].getTBBinom().getParameter();
		return this;
	}
	
	public boolean set(ModelStructure[][][] models, int t, int i) {
		// either all or no model exists (i.e. the label is not there for the condition ,e.g. no4sU)
		int present = 0;
		for (int s=0; s<buf.length; s++) {
			if (models[s][t][i]!=null) {
				System.arraycopy(models[s][t][i].getTBBinom().getParameter(), 0, buf[s], 0, buf[s].length);
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
			adaptBuf(buf[i],param);
			re+=tbbinom.logLik(data[i], buf[i]);
		}
		return re;
	}

	protected abstract void adaptBuf(double[] buf, double[] param);


	


}
