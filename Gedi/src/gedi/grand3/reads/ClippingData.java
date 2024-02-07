package gedi.grand3.reads;

import java.io.File;
import java.io.IOException;
import java.util.logging.Logger;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.correlation.Covariance;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import gedi.util.ArrayUtils;
import gedi.util.functions.EI;
import jdistlib.ChiSquare;

public class ClippingData {

	private int clip5p1;
	private int clip5p2;
	private int clip3p1;
	private int clip3p2;
	
	public ClippingData(int clip5p1, int clip5p2, int clip3p1, int clip3p2) {
		this.clip5p1 = clip5p1;
		this.clip5p2 = clip5p2;
		this.clip3p1 = clip3p1;
		this.clip3p2 = clip3p2;
	}
	
	
	public int getClip5p1() {
		return clip5p1;
	}
	
	public int getClip5p2() {
		return clip5p2;
	}
	
	public int getClip3p1() {
		return clip3p1;
	}
	
	public int getClip3p2() {
		return clip3p2;
	}
	
	@Override
	public String toString() {
		return "Read 1: "+clip5p1+","+clip3p1+" Read 2: "+clip5p2+","+clip3p2;
	}

	public static ClippingData fromFile(File file) throws IOException {
		int[] clip = EI.lines(file).splitField('\t',1).parseInt().skip(4).toIntArray();
		return new ClippingData(clip[0], clip[1],clip[2],clip[3]);
	}
	
	
	/**
	 * all position up to the one where the (bonferroni corrected) chisq pvalue of the mahalanobis distance is <0.01
	 * the mahalanobis distance is computed based on the mean and cov of the 20%-80% positions (running differences)
	 * 
	 * 
	 * @param m
	 * @return
	 */
	public static int[] inferClipping(double[][] m, Logger log) {
		if (m==null || m.length==0) return new int[] {0,0};

		for (int i=1; i<m.length; i++)
			m[i-1] = ArrayUtils.subtract(m[i-1],m[i]);
		m = ArrayUtils.redimPreserve(m, m.length-1);
			
		int s = (int) (ArrayUtils.nrows(m)*0.2);
		int e = (int) (ArrayUtils.nrows(m)*0.8);

		double[] mean = new double[ArrayUtils.ncols(m)];
		for (int d=0; d<mean.length; d++) 
			mean[d] = new Mean().evaluate(ArrayUtils.col(m,d),s,e-s);

		double[][] centered = new double[m.length][m[0].length];
		for (int i=0; i<centered.length; i++)
			for (int d=0; d<centered[0].length; d++)
				centered[i][d]=m[i][d]-mean[d];
		
		RealMatrix icov;
		try {
			icov = MatrixUtils.inverse(new Covariance(ArrayUtils.rows(m,s,e)).getCovarianceMatrix());
		} catch (SingularMatrixException ex) {
			log.severe("Cannot infer clipping parameters, matrix is singular!");
			log.severe(new Covariance(ArrayUtils.rows(m,s,e)).getCovarianceMatrix().toString());
			return new int[] {0,0};
		}
		double[] sqdist = new double[ArrayUtils.nrows(m)];
		int dim = mean.length;
		for (int i=0; i<sqdist.length; i++) 
			for (int r=0; r<dim; r++)
				for (int c=0; c<dim; c++)
					sqdist[i] += centered[i][r]*icov.getEntry(r, c)*centered[i][c];
		
		double[] p = sqdist;
		for (int i=0; i<p.length; i++) {
			p[i] = ChiSquare.cumulative(sqdist[i], mean.length, false, false);
			p[i] = Math.min(p[i]*p.length, 1);
		}
		
		int max = Math.min(25, p.length);
		int max2 = Math.min(10, p.length);
		int[] re = {0,0};
		for (int i=0; i<max; i++)
			if (p[i]<0.01)
				re[0] = i+1;
		for (int i=p.length-1; i>=p.length-max2; i--)
			if (p[i]<0.01)
				re[1] = p.length-i;
		if (re[0]==max || re[1]==max2) log.warning("Extreme clipping parameters, check!");

		return re;
		
	}


	public boolean isNoClipping() {
		return clip3p1==0 && clip5p1==0 && clip3p2==0 && clip5p2==0;
	}

	
}
