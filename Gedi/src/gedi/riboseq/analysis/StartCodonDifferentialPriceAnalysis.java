package gedi.riboseq.analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.startup.PriceStartUp;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.collections.PositionIterator;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.distributions.LfcDistribution;
import gedi.util.mutable.MutablePair;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;
import gedi.util.r.RRunner;

public class StartCodonDifferentialPriceAnalysis implements PriceAnalysis<Void> {

	/**
	 * For auto inclusion in {@link PriceStartUp}
	 */
	public static final String name = "differential-start-codon";
	public static final Function<GediParameterSet,GediParameter[]> params = set->new GediParameter[] {
		new GediParameter<Integer>(set, name+"-len", "Length after the start codon to consider", false, new IntParameterType(),100),
		new GediParameter<String>(set,name+"-contrasts", "Contrasts to compute changes around start codon", false, new StringParameterType(),"")
		
	};

	private int len;
	private LinkedHashMap<String, MutablePair<ToDoubleFunction<NumericArray>, ToDoubleFunction<NumericArray>>> contrasts;
	
	
	public StartCodonDifferentialPriceAnalysis(String[] conditions, PriceParameterSet param) throws IOException {
		this.len = (Integer)param.get(name+"-len").get();
		
		
		String contr = (String) param.get(name+"-contrasts").get();
		if (contr.length()==0) contr=null;
		
		HashMap<String, Integer> condIndex = ArrayUtils.createIndexMap(conditions);
		
		LinkedHashMap<String, MutablePair<ToDoubleFunction<NumericArray>, ToDoubleFunction<NumericArray>>> contrasts = new LinkedHashMap<>();
		
		if (contr==null && conditions.length==2)
			contrasts.put(conditions[0]+"/"+conditions[1], new MutablePair<>(a->a.getDouble(0), a->a.getDouble(1)));
		else if (contr!=null)
			for (String[] a : EI.lines(contr).skip(1).map(s->StringUtils.split(s, '\t')).loop()) 
				contrasts.put(a[0], new MutablePair<ToDoubleFunction<NumericArray>, ToDoubleFunction<NumericArray>>(createContrastPart(a[1],condIndex),createContrastPart(a[2],condIndex)));
		
		this.contrasts = contrasts.size()>0?contrasts:null;
	}
	
	private static ToDoubleFunction<NumericArray> createContrastPart(String d, HashMap<String, Integer> condIndex) {
		int[] ind = EI.split(d, '+').mapToDouble(condIndex::get).toIntArray();
		if (ind.length==0)
			return a->a.getDouble(ind[0]);
		else
			return a->{
				double re = 0;
				for (int i : ind)
					re+=a.getDouble(i);
				return re;
			};
	}

	@Override
	public void process(MajorIsoform data, LineWriter out, Void ctx) {
		//data.getTranscript().getData().getGeneId()
		if (contrasts!=null) {
			NumericArray sum = data.getSum(0);
			
			for (String co : contrasts.keySet()) {
				MutablePair<ToDoubleFunction<NumericArray>, ToDoubleFunction<NumericArray>> ctr = contrasts.get(co);
				double shift = LfcDistribution.mean(ctr.Item1.applyAsDouble(sum)+1,ctr.Item2.applyAsDouble(sum)+1);
				
				PositionIterator<NumericArray> it = data.iterateAminoAcids(0);
				while (it.hasNext()) {
					int p = it.nextInt();
					if (p>=len) 
						break;
					NumericArray a = it.getData();
					
					double mean;
					double var;
					if (a==null) {
						mean = LfcDistribution.mean(1,1);
						var = LfcDistribution.var(1,1);
					} else {
						mean = LfcDistribution.mean(ctr.Item1.applyAsDouble(a)+1,ctr.Item2.applyAsDouble(a)+1);
						var = LfcDistribution.var(ctr.Item1.applyAsDouble(a)+1,ctr.Item2.applyAsDouble(a)+1);
					}
					out.writef2("%s\t%s\t%d\t%.3f\t%.3f\n", data.getTranscript().getData().getGeneId(), co,p,mean-shift,var);
				}
			}
		}
	}

	@Override
	public void header(LineWriter out) throws IOException {
		if (contrasts!=null)
			out.writeLine("Gene\tContrast\tPosition\tMean\tVariance");
	}

	@Override
	public Void createContext() {
		return null;
	}

	@Override
	public void reduce(List<Void> ctx, LineWriter out) {
	}

	@Override
	public void plot(String data, String prefix) throws IOException {
			
		RRunner r = new RRunner(FileUtils.getFullNameWithoutExtension(data)+".R");
		r.set("prefix",prefix);
		r.set("output",FileUtils.getFullNameWithoutExtension(data)+".pdf");
		r.set("input",data);
		r.addSource(getClass().getResourceAsStream("/resources/R/differentialstart.R"));
		r.run(true);
	}

	
	
}
