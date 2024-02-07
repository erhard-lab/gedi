package gedi.core.data.annotation;


import java.io.IOException;

import gedi.app.Config;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class ScoreNameAnnotation implements BinarySerializable, ScoreProvider, NameProvider {
	
	private String name;
	private double score;
	
	public ScoreNameAnnotation() {
	}
	
	public ScoreNameAnnotation(String name, double score) {
		this.name = name;
		this.score = score;
	}
	

	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putString(name);
		out.putDouble(score);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		name = in.getString();
		score = in.getDouble();
	}
	
	@Override
	public String toString() {
		return String.format("%s("+Config.getInstance().getRealFormat()+")",name,score);
	}

	@Override
	public double getScore() {
		return score;
	}
	
	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public void setScore(double score) {
		this.score = score;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		long temp;
		temp = Double.doubleToLongBits(score);
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
		ScoreNameAnnotation other = (ScoreNameAnnotation) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (Double.doubleToLongBits(score) != Double
				.doubleToLongBits(other.score))
			return false;
		return true;
	}
	
	

}
