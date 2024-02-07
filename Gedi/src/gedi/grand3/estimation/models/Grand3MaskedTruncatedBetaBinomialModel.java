package gedi.grand3.estimation.models;

import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class Grand3MaskedTruncatedBetaBinomialModel implements MaximumLikelihoodModel<KNMatrixElement[]> {
	
	private double perr;
	private int pindex=0;
	private int shapeindex = 1;
	
	public Grand3MaskedTruncatedBetaBinomialModel(double perr) {
		this.perr = perr;
	}
	
	public Grand3MaskedTruncatedBetaBinomialModel(double perr, int pindex, int shapeindex) {
		this.perr = perr;
		this.pindex = pindex;
		this.shapeindex = shapeindex;
	}

	public double logLik(KNMatrixElement[] data, double[] param) {
		if (data.length==0) return 0;

		double p = param[pindex];
		double shape = param[shapeindex];
		
		double re = 0;
		int from = 0;
		for (int to=0; to<=data.length; to++) {
			
			if (to==data.length || data[to].n!=data[from].n) {
				int n=data[from].n;
				double norm = 0;
				double sum = 0;
				
				for (int i=from; i<to; i++) { 
					int k = data[i].k;
					re+=Grand3Model.ldtbbinom(k, n, perr,p,shape)*data[i].count;
					sum+=data[i].count;
					norm=norm==0?Grand3Model.ldtbbinom(k, n, perr,p,shape):Grand3Model.lse(norm, Grand3Model.ldtbbinom(k, n, perr,p,shape));
					// "from" is the first non-masked k, all up to n have to be added to norm!
					int nn = i==to-1?n:data[i+1].k-1;
					for (int kp=k+1; kp<=nn; kp++)
						norm=Grand3Model.lse(norm, Grand3Model.ldtbbinom(kp, n, perr,p,shape));
				}
				if (sum>0) 
					re-=norm*sum;
				
				from = to;
			}
		}
		return re;
	}

	
	public MaximumLikelihoodParametrization createPar() {
		return createPar(0.05,4.5);
	}

	public MaximumLikelihoodParametrization createPar(double pmconv, double shape) {
		return new MaximumLikelihoodParametrization(
				"maskedbinom",
				new String[] {"p.mconv","shape"},
				new double[] {pmconv,shape},
				new String[] {"%.4f","%.2f"},
				new double[] {0,-2},
				new double[] {1,5}
			);
	}
}