package gedi.slam;

import java.util.Arrays;
import java.util.logging.Logger;

import org.apache.commons.math3.special.Gamma;

import gedi.util.ArrayUtils;

public class VariationalBayesProportion {

	private double alpha;
	private double beta;
	
	private double p;
	private double q;
	
	private double thresVB = 1E-5;
	private double thresBinom = 1E-5;
	private int maxiter = 100;
	
//	private double[] binomest = new double[4];
	
	private double delta = Double.POSITIVE_INFINITY;
	
	/**
	 * 
	 * @param alpha Prior
	 * @param beta Prior
	 * @param p Mismatch prob
	 * @param q Conversion prob
	 */
	public VariationalBayesProportion(double alpha, double beta, double p, double q) {
		this.alpha = alpha;
		this.beta = beta;
		this.p = p;
		this.q = q;
	}
	
//	public boolean isPqConvergence() {
//		return delta<thresBinom;
//	}
//
//
//	
//	public void useEstimates() {
//		delta = Math.max(Math.abs(p-getEstimatedP()), Math.abs(q-getEstimatedQ()));
//		p = getEstimatedP();
//		q = getEstimatedQ();
//		Arrays.fill(binomest, 0);
//	}


	public double[] infer(int[] d, int[] n, double[] m, Logger log) {
		double[][] r = new double[2][d.length];
		
		// compute A and B
		double[] A = new double[d.length];
		double[] B = new double[d.length];
		for (int i=0; i<A.length; i++) {
			A[i] = d[i]*Math.log(q)+(n[i]-d[i])*Math.log(1-q);
			B[i] = d[i]*Math.log(p)+(n[i]-d[i])*Math.log(1-p);
		}
		
		double[] re = {0,0};
		double[] last = {Double.NEGATIVE_INFINITY,Double.NEGATIVE_INFINITY};
		int iter = 0;
		while (iter<maxiter && (Math.abs(re[0]-last[0])/re[0]>thresVB || Math.abs(re[1]-last[1])/re[1]>thresVB)) {
			double s1 = ArrayUtils.sum(r[0]);
			double s2 = ArrayUtils.sum(r[1]);
			double nor = Gamma.digamma(alpha+beta+s1+s2);
			
			double C = Gamma.digamma(alpha+s1)-nor;
			double D = Gamma.digamma(beta+s2)-nor;
			
			last[0] = re[0];
			last[1] = re[1];
			re[0] = alpha;
			re[1] = beta;
			
			for (int i=0; i<r[0].length; i++) {
				r[0][i] = Math.exp(A[i]+C);
				r[1][i] = Math.exp(B[i]+D);
				double s = r[0][i]+r[1][i];
				r[0][i] = r[0][i]/s;
				r[1][i] = r[1][i]/s;
				re[0]+=r[0][i]*m[i];
				re[1]+=r[1][i]*m[i];
			}
		}
		
		if (iter==maxiter) {
			log.warning("Variational Bayes did not converge within "+maxiter+" iterations (delta="+Math.max(Math.abs(re[0]-last[0])/re[0],Math.abs(re[1]-last[1])/re[1])+")!");
		}
		
//		synchronized (binomest) {
//			for (int i=0; i<d.length; i++) {
//				binomest[0]+=r[0][i]*d[i]*m[i];
//				binomest[1]+=r[0][i]*n[i]*m[i];
//				binomest[2]+=r[1][i]*d[i]*m[i];
//				binomest[3]+=r[1][i]*n[i]*m[i];
//			}
//		}
		
		return re;
	}
	
//	public double getEstimatedQ() {
//		return binomest[0]/binomest[1];
//	}
//	public double getEstimatedP() {
//		return binomest[2]/binomest[3];
//	}
//	public double getEstimatedTotal() {
//		return (binomest[0]+binomest[2])/(binomest[1]+binomest[3]);
//	}


	


	
	
}
