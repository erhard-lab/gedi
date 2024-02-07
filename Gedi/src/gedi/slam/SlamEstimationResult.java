package gedi.slam;

public class SlamEstimationResult {

	
	private double lower = 0;
	private double upper = 1;
	
	private double mean = Double.NaN;
	private double map = Double.NaN;
	
	private double alpha = Double.NaN;
	private double beta = Double.NaN;
	
	
	public SlamEstimationResult() {
	}
	
	public SlamEstimationResult(double lower, double mean, double map, double upper) {
		this.lower = lower;
		this.upper = upper;
		this.map = map;
		this.mean = mean;
	}
	
	public SlamEstimationResult(double lower, double mean, double map, double upper, double alpha, double beta) {
		this.lower = lower;
		this.upper = upper;
		this.map = map;
		this.mean = mean;
		this.alpha = alpha;
		this.beta = beta;
	}

	public double getLower() {
		return lower;
	}

	public double getUpper() {
		return upper;
	}

	public double getMean() {
		return mean;
	}

	public double getMap() {
		return map;
	}

	public double getAlpha() {
		return alpha;
	}

	public double getBeta() {
		return beta;
	}
	
	
	
}
