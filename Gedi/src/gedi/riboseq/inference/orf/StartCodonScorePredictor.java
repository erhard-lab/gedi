package gedi.riboseq.inference.orf;

import java.io.IOException;

import gedi.riboseq.inference.orf.StartCodonTraining.RangeAndStartClassifier;
import gedi.riboseq.inference.orf.StartCodonTraining.StartCodonPredictionOrf;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class StartCodonScorePredictor implements BinarySerializable {

	private RangeAndStartClassifier model;
	private int upstream;
	private int downstream;
	
	public StartCodonScorePredictor() {
	}
	
	public StartCodonScorePredictor(RangeAndStartClassifier model, int upstream, int downstream) {
		this.model = model;
		this.upstream = upstream;
		this.downstream = downstream;
	}
	
	
	public void predict(PriceOrf orf) {
		StartCodonPredictionOrf o = new StartCodonPredictionOrf(orf,0);
		double[] range = new double[o.length()];
		double[] start = new double[o.length()];
		double[] x = null;
		double[] buff = new double[2];
		for (int p=0; p<o.length(); p++) {
			x = StartCodonTraining.getX(o, p, upstream, downstream, x);
			range[p] = model.predictRange(x, buff);
			start[p] = model.predictStart(x, buff);
		}
		orf.startScores = start;
		orf.startRangeScores = range;
	}


	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putCInt(upstream);
		out.putCInt(downstream);
		model.serialize(out);
	}


	@Override
	public void deserialize(BinaryReader in) throws IOException {
		upstream = in.getCInt();
		downstream = in.getCInt();
		model = new RangeAndStartClassifier();
		model.deserialize(in);
	}
	
}
