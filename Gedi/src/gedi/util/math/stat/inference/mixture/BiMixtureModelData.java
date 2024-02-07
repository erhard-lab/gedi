package gedi.util.math.stat.inference.mixture;

public class BiMixtureModelData {

	double ll0;
	double ll1;
	double count;
	
	
	public BiMixtureModelData(double ll0, double ll1, double count) {
		super();
		this.ll0 = ll0;
		this.ll1 = ll1;
		this.count = count;
	}

	public BiMixtureModelData set(double ll0, double ll1, double count) {
		this.ll0 = ll0;
		this.ll1 = ll1;
		this.count = count;
		return this;
	}
	

	public boolean isInt() {
		return (int)count==count;
	}


	@Override
	public String toString() {
		return "[" + ll0 + "," + ll1 + "]x"+count;
	}


	public double getCount() {
		return count;
	}
	
	
	
	
}
