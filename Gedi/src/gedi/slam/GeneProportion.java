package gedi.slam;

import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;

public class GeneProportion {

	private GeneData data;
	private SlamEstimationResult[] est;
	public GeneProportion(GeneData data, SlamEstimationResult[] est) {
		this.data = data;
		this.est = est;
	}

	public String getGene() {
		return data.getGene();
	}

	public GeneData getData() {
		return data;
	}
	
	
	public NumericArray getConversions() {
		return data.getTotalConversions();
	}

	public NumericArray getCoverage() {
		return data.getTotalCoverage();
	}

	public NumericArray getDoubleHits() {
		return data.getTotalDoubleHits();
	}

	public NumericArray getDoubleHitCoverage() {
		return data.getTotalDoubleHitCoverage();
	}

	public NumericArray getReadCount() {
		return data.getReadCount();
	}
	
	public SlamEstimationResult getEstimated(int cond) {
		return est[cond];
	}
	
	@Override
	public String toString() {
		return StringUtils.toString(est);
	}
	
}
