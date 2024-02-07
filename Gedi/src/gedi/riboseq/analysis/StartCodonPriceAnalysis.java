package gedi.riboseq.analysis;

import java.io.IOException;
import java.util.function.Function;

import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.startup.PriceStartUp;
import gedi.util.FileUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.collections.PositionIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.r.RRunner;

public class StartCodonPriceAnalysis extends MetaGenePriceAnalysis<Integer> {

	/**
	 * For auto inclusion in {@link PriceStartUp}
	 */
	public static final String name = "start-codon";
	public static final Function<GediParameterSet,GediParameter[]> params = set->new GediParameter[] {
		new GediParameter<Integer>(set, name+"-len", "Length after the start codon to consider", false, new IntParameterType(),50)	
	};

	private int len;
	
	
	public StartCodonPriceAnalysis(String[] conditions, PriceParameterSet param) {
		super(conditions,"Position");
		this.len = (Integer)param.get(name+"-len").get();
	}

	@Override
	public void process(MajorIsoform data, LineWriter out, NormalizingCounter<Integer>[] ctx) {
		
		NumericArray sum = data.getSum(0);
		
		PositionIterator<NumericArray> it = data.iterateAminoAcids(0);
		while (it.hasNext()) {
			int p = it.nextInt();
			if (p>=len) 
				break;
			NumericArray a = it.getData();
			for (int i=0; i<conditions.length; i++){
				ctx[i].count(p, a==null||sum.getDouble(i)==0?0:(a.getDouble(i)/sum.getDouble(i)*data.getAminoAcidLength(0)));
			}
		}
	}


	@Override
	public void plot(String data, String prefix) throws IOException {
		RRunner r = new RRunner(FileUtils.getFullNameWithoutExtension(data)+".R");
		r.set("prefix",prefix);
		r.set("output",FileUtils.getFullNameWithoutExtension(data)+".png");
		r.set("input",data);
		r.addSource(getClass().getResourceAsStream("/resources/R/around.R"));
		r.run(true);
	}
	
}
