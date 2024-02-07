package gedi.grand3.estimation;

import gedi.grand3.estimation.models.Grand3Model;

public class OldNewLogLikCache {

	
	private static final int PRE = 100;
	private double[][] tabOld;
	private double[][] tabNew;
	private Grand3Model model;
	private double[] param;
	
	
	public OldNewLogLikCache(Grand3Model model, double[] param) {
		this(model,param, PRE);
	}
	public OldNewLogLikCache(Grand3Model model, double[] param ,int pre) {
		this.model = model;
		this.param = param;
		
		tabOld = new double[pre][pre];
		tabNew = new double[pre][pre];
		for (int k=0; k<pre; k++)
			for (int n=0; n<pre; n++) {
				tabOld[k][n] = model.logLikOld(k, n, param); 
				tabNew[k][n] = model.logLikNew(k, n, param); 
			}
	}
	
	public double getNewLogLik(int k, int n) {
		if (n>=PRE)
			return model.logLikNew(k, n, param);
		return tabNew[k][n];
	}

	public double getOldLogLik(int k, int n) {
		if (n>=PRE)
			return model.logLikOld(k, n, param);
		return tabOld[k][n];
	}

	
}
