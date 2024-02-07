package gedi.util.math.stat.inference.mixture;

import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class BiMixtureModelResult implements BinarySerializable{
	private double lower;
	private double map;
	private double upper;
	private double alpha;
	private double beta;
	public BiMixtureModelResult() {}
	public BiMixtureModelResult(double lower, double map, double upper, double alpha, double beta) {
		if (Double.isNaN(lower) || Double.isNaN(map) || Double.isNaN(upper) || lower<0 || map<0 || upper<lower || map>1 || upper>1)
			throw new RuntimeException("Invalid parameters!");
		this.lower = lower;
		this.map = map;
		this.upper = upper;
		this.alpha = alpha;
		this.beta = beta;
	}
	public double getLower() {
		return lower;
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
	@Override
	public String toString() {
		return "lower=" + lower + ", map=" + map + ", upper=" + upper + ", alpha=" + alpha + ", beta="
				+ beta ;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(alpha);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(beta);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(lower);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(map);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(upper);
		result = prime * result + (int) (temp ^ (temp >>> 32));
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
		if (Double.doubleToLongBits(alpha) != Double.doubleToLongBits(other.alpha))
			return false;
		if (Double.doubleToLongBits(beta) != Double.doubleToLongBits(other.beta))
			return false;
		if (Double.doubleToLongBits(lower) != Double.doubleToLongBits(other.lower))
			return false;
		if (Double.doubleToLongBits(map) != Double.doubleToLongBits(other.map))
			return false;
		if (Double.doubleToLongBits(upper) != Double.doubleToLongBits(other.upper))
			return false;
		return true;
	}
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putFloat(lower);
		out.putFloat(map);
		if (Double.isNaN(alpha)) {
			out.putFloat(-upper);
		}
		else {
			out.putFloat(upper);
			out.putFloat(alpha);
			out.putFloat(beta);
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
		}
		else {
			alpha = in.getFloat();
			beta = in.getFloat();
		}
	}
	
	
}
