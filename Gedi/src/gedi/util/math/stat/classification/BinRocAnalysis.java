package gedi.util.math.stat.classification;

import gedi.util.ArrayUtils;

import java.util.Arrays;

public class BinRocAnalysis extends RocAnalysis {

	private int[] t;
	private int[] f;
	
	private int[] tp;
	private int[] fp;
	private double[] bins;
	
	private boolean computed = false;
	private boolean reverseFlag;
	
	/**
	 * Wheter or not to invert meaning of score
	 * @param reverseFlag
	 */
	public BinRocAnalysis(boolean reverseFlag, double...bins) {
		this.reverseFlag = reverseFlag;
		t = new int[bins.length+1];
		f = new int[bins.length+1];
		this.bins = bins;
		Arrays.sort(bins);
	}

	@Override
	/**
	 * The larger the score, the better!
	 */
	public void addScore(double score, boolean isTrue) {
		int bin = Arrays.binarySearch(bins, score);
		if (bin<0) bin = -bin-1;
		bin = Math.min(bin, t.length-1);
		if (isTrue)
			t[bin]++;
		else
			f[bin]++;
		computed=false;
	}
	
	@Override
	public void switchDirection() {
		reverseFlag=!reverseFlag;
		computed = false;
	}
	
	@Override
	public int getNumPoints() {
		return t.length-1;
	}
	
	public int getNumPositives() {
		compute();
		return tp[reverseFlag?tp.length-1:0];
	}
	
	public int getNumNegatives() {
		compute();
		return fp[reverseFlag?fp.length-1:0];
	}
	
	protected void compute() {
		if (computed)return;
		computed = true;
		tp=ArrayUtils.cumSum(t, reverseFlag?1:-1);
		fp=ArrayUtils.cumSum(f, reverseFlag?1:-1);
	}

	@Override
	public double getCutoff(int index) {
		return bins[index];
	}

	@Override
	public int getFn(int index) {
		compute();
		return getNumPositives()-tp[index];
	}

	@Override
	public int getFp(int index) {
		compute();
		return fp[index];
	}

	@Override
	public int getTn(int index) {
		compute();
		return getNumNegatives()-fp[index]; 
	}

	@Override
	public int getTp(int index) {
		compute();
		return tp[index];
	}

	
}
