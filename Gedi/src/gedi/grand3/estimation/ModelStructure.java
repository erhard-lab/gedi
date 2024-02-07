package gedi.grand3.estimation;

import java.io.IOException;
import java.util.HashMap;

import gedi.grand3.experiment.ExperimentalDesign;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class ModelStructure implements BinarySerializable, Comparable<ModelStructure> {

	private int s;
	private int t;
	private int i;
	private MaximumLikelihoodParametrization binom;
	private MaximumLikelihoodParametrization tbbinom;
	
	/**
	 * For deserialization, provide empty models
	 * @param binom
	 * @param tbbinom
	 */
	public ModelStructure(MaximumLikelihoodParametrization binom,
			MaximumLikelihoodParametrization tbbinom) {
		this.binom = binom;
		this.tbbinom = tbbinom;
	}
	
	public ModelStructure(int s, int t, int i, MaximumLikelihoodParametrization binom,
			MaximumLikelihoodParametrization tbbinom) {
		this.s = s;
		this.t = t;
		this.i = i;
		this.binom = binom;
		this.tbbinom = tbbinom;
	}


	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putCInt(s);
		out.putCInt(t);
		out.putCInt(i);
		binom.serialize(out);
		tbbinom.serialize(out);
	}


	@Override
	public void deserialize(BinaryReader in) throws IOException {
		s=in.getCInt();
		t=in.getCInt();
		i=in.getCInt();
		binom.deserialize(in);
		tbbinom.deserialize(in);
	}
	
	public int getSubread() {
		return s;
	}
	public int getType() {
		return t;
	}
	public int getSample() {
		return i;
	}
	
	public void writeTableHeader(LineWriter out, boolean ci) throws IOException {
		if (ci)
			out.writef("Estimator\tCondition\tSubread\tLabel\tLower prior p.err\tUpper prior p.err\t%s\t%s\t%s\tBinom log likelihood\t%s\t%s\t%s\tTB-Binom log likelihood\n",
					EI.wrap(binom.getParamNames()).map(c->"Binom "+c).concat("\t"),
					EI.wrap(binom.getParamNames()).map(c->"Lower Binom "+c).concat("\t"),
					EI.wrap(binom.getParamNames()).map(c->"Upper Binom "+c).concat("\t"),
					EI.wrap(tbbinom.getParamNames()).map(c->"TB-Binom "+c).concat("\t"),
					EI.wrap(tbbinom.getParamNames()).map(c->"Lower TB-Binom "+c).concat("\t"),
					EI.wrap(tbbinom.getParamNames()).map(c->"Upper TB-Binom "+c).concat("\t"));
		else
			out.writef("Estimator\tCondition\tSubread\tLabel\tLower prior p.err\tUpper prior p.err\t%s\tBinom log likelihood\t%s\tTB-Binom log likelihood\n",
					EI.wrap(binom.getParamNames()).map(c->"Binom "+c).concat("\t"),
					EI.wrap(tbbinom.getParamNames()).map(c->"TB-Binom "+c).concat("\t"));
	}

	public void writeTable(LineWriter out, String estimator, boolean ci, ExperimentalDesign design, String[] subreads, double[][][][] pre_perr) throws IOException {
		out.writef("%s\t%s\t%s\t%s\t%.4g\t%.4g\t",estimator,design.getSampleNameForSampleIndex(i),subreads[s],design.getTypes()[t].toString(),pre_perr[s][t][i][0],pre_perr[s][t][i][1]);
		out.write(binom.toTabString());
		if (ci){
			out.write("\t");
			out.write(binom.toTabStringLower());
			out.write("\t");
			out.write(binom.toTabStringUpper());
		}
		out.writef("\t%.2f\t", binom.getLogLik());
		out.write(tbbinom.toTabString());
		if (ci){
			out.write("\t");
			out.write(tbbinom.toTabStringLower());
			out.write("\t");
			out.write(tbbinom.toTabStringUpper());
		}
		out.writef("\t%.2f\n", tbbinom.getLogLik());
	}
	public void parse(String[] a, HeaderLine h, HashMap<String, Integer> sampleIndex, HashMap<String, Integer> subreadIndex, HashMap<String, Integer> typeIndex) {
		i=sampleIndex.get(a[h.get("Condition")]);
		s=subreadIndex.get(a[h.get("Subread")]);
		t=typeIndex.get(a[h.get("Label")]);
		
		binom.parse(a,h,"Binom");
		tbbinom.parse(a,h,"TB-Binom");
	}
	
	
	
	public MaximumLikelihoodParametrization getBinom() {
		return binom;
	}
	public MaximumLikelihoodParametrization getTBBinom() {
		return tbbinom;
	}
	@Override
	public int compareTo(ModelStructure o) {
		int re = Integer.compare(t, o.t);
		if (re!=0) return re;
		re = Integer.compare(i, o.i);
		if (re!=0) return re;
		return Integer.compare(s, o.s);
	}
	
	
	public static int compareNoS(ModelStructure a, ModelStructure b) {
		int re = Integer.compare(a.t, b.t);
		if (re!=0) return re;
		return Integer.compare(a.i, b.i);
	}

	
	
	
	
}
