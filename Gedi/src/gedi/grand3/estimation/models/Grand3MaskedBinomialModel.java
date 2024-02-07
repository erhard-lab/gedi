package gedi.grand3.estimation.models;

import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class Grand3MaskedBinomialModel implements MaximumLikelihoodModel<KNMatrixElement[]> {
	
	private int pindex = 0;
	
	public Grand3MaskedBinomialModel() {
	}
	public Grand3MaskedBinomialModel(int pindex) {
		this.pindex = pindex;
	}


	public double logLik(KNMatrixElement[] data, double[] param) {
		if (data.length==0) return 0;
		
		double p = param[pindex];
		
		double re = 0;
		int from = 0;
		for (int to=0; to<=data.length; to++) {
			
			if (to==data.length || data[to].n!=data[from].n) {
				int n=data[from].n;
				double norm = 0;
				double sum = 0;
				
				for (int i=from; i<to; i++) { 
					int k = data[i].k;
					re+=Grand3Model.lbinom(k, n, p)*data[i].count;
					sum+=data[i].count;
					norm=norm==0?Grand3Model.lbinom(k, n, p):Grand3Model.lse(norm, Grand3Model.lbinom(k, n, p));
					// "from" is the first non-masked k, all up to n have to be added to norm!
					int nn = i==to-1?n:data[i+1].k-1;
					for (int kp=k+1; kp<=nn; kp++)
						norm=Grand3Model.lse(norm, Grand3Model.lbinom(kp, n, p));
				}
				if (sum>0) 
					re-=norm*sum;
				from = to;
			}
		}
		return re;
	}

	
	public MaximumLikelihoodParametrization createPar() {
		return createPar(0.05);
	}

	public MaximumLikelihoodParametrization createPar(double pconv) {
		return new MaximumLikelihoodParametrization(
				"maskedbinom",
				new String[] {"p.conv"},
				new double[] {pconv},
				new String[] {"%.4f"},
				new double[] {0},
				new double[] {1}
			);
	}


}