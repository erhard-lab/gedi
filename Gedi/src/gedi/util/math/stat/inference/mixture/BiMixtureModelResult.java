package gedi.util.math.stat.inference.mixture;

import java.io.IOException;
import java.util.Arrays;
import java.util.Objects;

import gedi.util.ArrayUtils;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class BiMixtureModelResult implements BinarySerializable{
	private double lower;
	private double map;
	private double upper;
	private double alpha;
	private double beta;
	private double integral;
	private double[] discrete0;
	
	public BiMixtureModelResult() {}
	public BiMixtureModelResult(double lower, double map, double upper, double alpha, double beta, double integral, double[] discrete0) {
		if (Double.isNaN(lower) || Double.isNaN(map) || Double.isNaN(upper) || lower<0 || map<0 || upper<lower || map>1 || upper>1)
			throw new RuntimeException("Invalid parameters!");
		this.lower = lower;
		this.map = map;
		this.upper = upper;
		this.alpha = alpha;
		this.beta = beta;
		this.integral = integral;
		this.discrete0 = discrete0;
	}
	public double getLower() {
		return lower;
	}
	
//	public double getBetaMixNtr() {
//		return mixBeta==null?Double.NaN:mixBeta[0];
//	}
//	public double getBetaMixMix() {
//		return mixBeta==null?Double.NaN:mixBeta[1];
//	}
//	public double getBetaMixAlpha1() {
//		return mixBeta==null?Double.NaN:mixBeta[2];
//	}
//	public double getBetaMixBeta1() {
//		return mixBeta==null?Double.NaN:mixBeta[3];
//	}
//	public double getBetaMixAlpha2() {
//		return mixBeta==null?Double.NaN:mixBeta[4];
//	}
//	public double getBetaMixBeta2() {
//		return mixBeta==null?Double.NaN:mixBeta[5];
//	}
//	public double getBetaMixSSE() {
//		return mixBeta==null?Double.NaN:mixBeta[6];
//	}
//	public double getBetaMixIntegral() {
//		return mixBeta==null?Double.NaN:mixBeta[7];
//	}
	
	public double[] getDiscrete0() {
		return discrete0;
	}
	public double getMap() {
		return map;
	}
	public double getUpper() {
		return upper;
	}
	public double getAlpha() {
		return alpha;
	}
	public double getBeta() {
		return beta;
	}
	public double getIntegral() {
		return integral;
	}
	@Override
	public String toString() {
		return "lower=" + lower + ", map=" + map + ", upper=" + upper + ", alpha=" + alpha + ", beta="
				+ beta + ", integral="+integral+(discrete0!=null?" discrete0="+Arrays.toString(discrete0):"");
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(discrete0);
		result = prime * result + Objects.hash(alpha, beta, integral, lower, map, upper);
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		BiMixtureModelResult other = (BiMixtureModelResult) obj;
		return Double.doubleToLongBits(alpha) == Double.doubleToLongBits(other.alpha)
				&& Double.doubleToLongBits(beta) == Double.doubleToLongBits(other.beta)
				&& Double.doubleToLongBits(integral) == Double.doubleToLongBits(other.integral)
				&& Double.doubleToLongBits(lower) == Double.doubleToLongBits(other.lower)
				&& Double.doubleToLongBits(map) == Double.doubleToLongBits(other.map)
				&& Arrays.equals(discrete0, other.discrete0)
				&& Double.doubleToLongBits(upper) == Double.doubleToLongBits(other.upper);
	}
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putFloat(discrete0==null?lower:-lower);
		out.putFloat(map);
		if (Double.isNaN(alpha)) {
			out.putFloat(-upper);
		}
		else {
			out.putFloat(upper);
			out.putFloat(alpha);
			out.putFloat(beta);
			out.putFloat(integral);
		}
		if (discrete0!=null) {
			out.putCInt(discrete0.length);
			for (int i=0; i<discrete0.length; i++)
				out.putFloat(discrete0[i]);
		}
	}
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		lower = in.getFloat();
		map = in.getFloat();
		upper = in.getFloat();
		if (1/upper<0) {
			// 1/ is necessary as upper might be 0!
			upper=-upper;
			alpha = Double.NaN;
			beta = Double.NaN;
			integral = Double.NaN;
		}
		else {
			alpha = in.getFloat();
			beta = in.getFloat();
			integral = in.getFloat();
		}
		if (1/lower<0) {
			lower = -lower;
			discrete0 = new double[in.getCInt()];
			for (int i=0; i<discrete0.length; i++)
				discrete0[i] = in.getFloat();
		} else
			discrete0 = null;
	}
	
	
}
