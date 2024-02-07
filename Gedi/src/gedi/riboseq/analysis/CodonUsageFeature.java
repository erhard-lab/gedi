package gedi.riboseq.analysis;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.ConcurrentModificationException;
import java.util.HashMap;
import java.util.Set;

import javax.script.ScriptException;

import gedi.core.genomic.Genomic;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.features.AbstractFeature;
import gedi.core.region.feature.output.Barplot;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.nashorn.JS;
import gedi.util.userInteraction.results.ResultProducer;
import org.openjdk.nashorn.api.scripting.ScriptObjectMirror;


public class CodonUsageFeature extends AbstractFeature<Void> {

	private HashMap<String,NumericArray[]> counters = new HashMap<String,NumericArray[]>();
	
	// temp
	private HashMap<String,NumericArray[]> buffers = new HashMap<String,NumericArray[]>();
	private NumericArray total;
	private NumericArray codonCount;
	
	private int decimals=2;
	private ArrayList<Barplot> plots = new ArrayList<Barplot>(); 
	
	private double minEachConditionThreshold = 0;
	private double sumConditionThreshold = 0;
	private int offsetFromStart = 15;
	private int offsetFromEnd = 15;
	
	private int beforeE = 0;
	private int afterA = 3;

	private double minThreshold = 1;
	private boolean maxDownsampling = false;
	private boolean normalizeExpression = false;
	private Genomic genomic;
	
	
	private ArrayList<CodonUsageOutputOption> outputs = new ArrayList<>();
	
	public CodonUsageFeature(String file) {
		minValues = maxValues = 0;
		minInputs = maxInputs = 0;
		setFile(file);
	}
	
	public void setFile(String path) {
		setId(path);
	}
	
	public void setGenomic(Genomic genomic) {
		this.genomic = genomic;
	}
	
	
	public void add(CodonUsageOutputOption option) {
		this.outputs.add(option);
	}
	
	public void addOutput(String label, String formula) throws ScriptException {
		add(new CodonUsageOutputOption(label, formula));
	}
	
	public void addRepresentationA() throws ScriptException {
		addOutput("A-site", "A/(d[0]+d[1]+d[2])*3");
		afterA = Math.max(afterA, 3);
	}
	
	public void addSticky() throws ScriptException {
		addOutput("Sticky-site", "u[2]");
		beforeE=9;
		afterA = Math.max(afterA, 3);
	}
	
	public void setMinEachConditionThreshold(double minEachConditionThreshold) {
		this.minEachConditionThreshold = minEachConditionThreshold;
	}

	public void setDecimals(int decimals) {
		this.decimals = decimals;
	}
	
	public int getDecimals() {
		return decimals;
	}
	
	public void setMaxDownsampling() {
		this.maxDownsampling = true;
	}
	
	public void setNormalize() {
		this.normalizeExpression = true;
	}
	
	public void add(Barplot plot) {
		plots.add(plot);
	}
	
	
	@Override
	public CodonUsageFeature addResultProducers(ArrayList<ResultProducer> re) {
		re.addAll(plots);
		return this;
	}
	
	
	
	@Override
	protected void accept_internal(Set<Void> values) {
		
		
		
		PriceOrf orf = (PriceOrf) referenceRegion.getData();
		int start = Math.max(orf.getPredictedStartAminoAcid(),0)+Math.max(beforeE+1, offsetFromStart);
		String seq = genomic.getSequence(referenceRegion).toString();
		seq = seq.substring(orf.getPredictedStartAminoAcid()*3);
		
		for (NumericArray[] c : buffers.values())
			for (NumericArray a : c)
				a.clear();
		
		
		int ncount = 3+beforeE+afterA;
		
		if (total==null)
			total = NumericArray.createMemory(orf.getNumConditions(), NumericArrayType.Double);
		else
			total.clear();
		
		if (codonCount==null)
			codonCount = NumericArray.createMemory(orf.getNumConditions(), NumericArrayType.Integer);
		else
			codonCount.clear();
		
//		int[] perm = EI.seq(0,act[0].length).toIntArray();
//		ArrayUtils.shuffleSlice(perm, start, perm.length);
		
		for (int i=start; i<orf.getOrfAaLength()-Math.max(offsetFromEnd,afterA+1); i++) {
			double downsamplingFactor = 1;
			if (maxDownsampling) {
				for (int c=0; c<orf.getNumConditions(); c++)
					downsamplingFactor = Math.max(downsamplingFactor,orf.getProfile(c, i));
				downsamplingFactor = 1/downsamplingFactor;
			}
			for (int c=0; c<orf.getNumConditions(); c++) {
				if (orf.getProfile(c, i)>=minThreshold ) {
					
					total.add(c,orf.getProfile(c, i)*downsamplingFactor);
					codonCount.add(c, 1);
					
					
					for (int off=0; off<ncount; off++) {
						int aapos = i+off-beforeE-1;
						String codon = seq.substring(aapos*3, aapos*3+3);
						NumericArray[] bi = buffers.computeIfAbsent(codon, x->createBuffer(orf.getNumConditions(),ncount));
						bi[off].add(c,orf.getProfile(c, i)*downsamplingFactor);
					}
					
				}
			}
		}
		
		
		if (NumericArrayFunction.Min.applyAsDouble(total)>=minEachConditionThreshold &&
				NumericArrayFunction.Sum.applyAsDouble(total)>=sumConditionThreshold) {

			if (normalizeExpression) {
				for (NumericArray[] arr : buffers.values())
					for (NumericArray a : arr)
						for (int c=0; c<a.length(); c++) {
							double norma = a.getDouble(c)/total.getDouble(c)*codonCount.getInt(c);
							if (Double.isNaN(norma))
								norma = 0;
							a.setDouble(c, norma);
						}
			}
			
	
			for (String codon : buffers.keySet()) {
				NumericArray[] cou = counters.computeIfAbsent(codon, x->createBuffer(orf.getNumConditions(), ncount));
				NumericArray[] buff = buffers.get(codon);
				for (int i=0; i<cou.length; i++)
					cou[i].add(buff[i]);
			}
		
		}
		
	}

	private NumericArray[] createBuffer(int eachlength, int count) {
		NumericArray[] re = new NumericArray[count];
		for (int i = 0; i < re.length; i++) {
			re[i] = NumericArray.createMemory(eachlength, NumericArrayType.Double);
		}
		return re;
	}

	@Override
	public GenomicRegionFeature<Void> copy() {
		CodonUsageFeature re = new CodonUsageFeature(getId());
		re.copyProperties(this);
		re.plots = plots;
		re.decimals = decimals;
		re.genomic = genomic;
		return re;
	}
	
	@Override
	public void produceResults(GenomicRegionFeature<Void>[] o) {
		if (o!=null) {
			
			for (NumericArray[] c : counters.values())
				for (NumericArray a : c)
					a.clear();
			
			for (GenomicRegionFeature<Void> a : o) {
				CodonUsageFeature x = (CodonUsageFeature)a;
				
				String[] codons = null;
				while (codons==null) {
					try {
						codons = x.counters.keySet().toArray(new String[0]);
					} catch (ConcurrentModificationException e) {
						codons = null;
					}
				}
				for (String c : codons) {
					NumericArray[] ot = x.counters.get(c);
					NumericArray[] mi = counters.computeIfAbsent(c, asdsd->createBuffer(ot[0].length(), ot.length));
					
					for (int i=0; i<mi.length; i++)
						mi[i].add(ot[i]);
				}
			}
			throw new RuntimeException("Not yet implemented!");
			
		}
		produceResult(!program.isRunning());
	}

	private void produceResult(boolean isFinal) {
		if (counters.size()==0) return;
		
		LineOrientedFile out = new LineOrientedFile(getId());
		try {
			out.startWriting();
			out.writef("Amino acid\tCodon\tStatistic");
			
			int l = counters.values().iterator().next()[0].length();
			if (program.getLabels()!=null && program.getLabels().length==l)
				for (int i=0; i<program.getLabels().length; i++) 
					out.writef("\t%s",program.getLabels()[i]);
			else {
				for (int i=0; i<l; i++) 
					out.writef("\t%d",i);
			}
			out.writeLine();
			
			
			HashMap<String,NumericArray[]> counters = this.counters;
			String[] codons = counters.keySet().toArray(new String[0]);
			
			Comparator<String> aaComp = (a,b)->SequenceUtils.translate(a).compareTo(SequenceUtils.translate(b));
			Arrays.sort(codons, aaComp.thenComparing((a,b)->a.compareTo(b)));
			
			for (String t : codons) {
				NumericArray[] a = counters.get(t);
				
				for (CodonUsageOutputOption opt : outputs) {
					out.writef("%s\t%s\t%s",SequenceUtils.translate(t),t,opt.label);
					for (int c=0; c<a[0].length(); c++)
						out.writef("\t%.5f",opt.compute(a,c, beforeE));
					out.writeLine();
				}
			}
			
			
			out.finishWriting();
			
			
			for (Barplot plot : plots) {
				plot.setName(getId());
				plot.plot(out,isFinal,new String[] {"Amino acid","Codon"},inputNames.length>1);
			}
			
			
		} catch (IOException e) {
			throw new RuntimeException("Cannot write output file!",e);
		}
	}
	
	
	public static class CodonUsageOutputOption {
		private String label;
		private ScriptObjectMirror p;
		// A/(d[0]+d[1]+d[2])*3

		public CodonUsageOutputOption(String label, String formula) throws ScriptException {
			this.label = label;
			
			StringBuilder code = new StringBuilder();
			if (formula.contains(";") || formula.contains("\n"))
				code.append("function(u,E,P,A,d) {\n").append(formula).append("}");
			else
				code.append("function(u,E,P,A,d) ").append(formula);
			p = new JS().execSource(code.toString());
		}

		
		public double compute(NumericArray[] counts, int condition, int beforeA) {
			int afterE = counts.length-beforeA-3;
			
			double[] u = beforeA>0?new double[beforeA]:null;
			double[] d = afterE>0?new double[afterE]:null;
			
			for (int i=0; i<beforeA; i++)
				u[beforeA-1-i] = counts[i].getDouble(condition);
			double E = counts[beforeA].getDouble(condition);
			double P = counts[beforeA+1].getDouble(condition);
			double A = counts[beforeA+2].getDouble(condition);
			for (int i=0; i<afterE; i++)
				d[i] = counts[beforeA+3+i].getDouble(condition);

			double re = (double) p.call(null, u,E,P,A,d);
			return re;
		}


	}
	
	

}
