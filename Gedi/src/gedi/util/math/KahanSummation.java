package gedi.util.math;

/**
 * Numerically stable summation
 * @author erhard
 *
 */
public class KahanSummation {

	private double sum = 0;
	private double c = 0;

	public void add(double s) {
		double y = s - c;
		double t = sum + y; 
		c = (t - sum) - y;
		sum = t;           
	}

	public double getSum() {
		return sum;
	}

}
