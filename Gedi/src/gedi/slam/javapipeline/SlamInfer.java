package gedi.slam.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.Function;
import java.util.logging.Level;
import java.util.regex.Pattern;

import gedi.core.genomic.Genomic;
import gedi.slam.GeneData;
import gedi.slam.OptimNumericalIntegrationProportion;
import gedi.slam.SlamEstimationResult;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.randomaccess.serialization.BinarySerializableSerializer;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediParameterSpec;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;
import gedi.util.userInteraction.progress.Progress;

public class SlamInfer extends GediProgram {

	
	public SlamInfer(SlamParameterSet params) {
		addInput(params.nthreads);
		addInput(params.dataFile);
		addInput(params.rateTable);
		addInput(params.lower);
		addInput(params.upper);
		addInput(params.biasnew);
		addInput(params.genomic);
		addInput(params.full);
		addInput(params.no4sUpattern);
		addInput(params.dataExtFile);
		addInput(params.sparseOutput);
		addInput(params.doubleOnly);
		addInput(params.prefix);
		
		addOutput(params.extTable);
	}
	
	
	@Override
	protected void initParameter(GediParameterSet param, GediParameterSpec inputSpec) {
		SlamParameterSet params = (SlamParameterSet)param;
		
		if (params.sparseOutput.get())
			addOutput(params.out10x);
		else
			addOutput(params.outTable);
	}
	
	
	public String execute(GediProgramContext context) throws IOException {
		
		int nthreads = getIntParameter(0);
		File data = getParameter(1);
		File rateFile = getParameter(2);
		double lower = getDoubleParameter(3);
		double upper = getDoubleParameter(4);
		double bias = getDoubleParameter(5);
		Genomic genomic = getParameter(6);
		boolean full = getParameter(7);
		String pat = getParameter(8);
		File ext = getParameter(9);
		boolean o10x = getParameter(10);
		boolean doubleonly = getParameter(11);
		String prefix = getParameter(12);
		
		ExtendedIterator<String[]> rit = EI.lines(rateFile.getPath()).map(s->StringUtils.split(s,'\t'));
		String[] cond=ArrayUtils.slice(rit.next(),1);
		
		
		Pattern no4sUPattern = Pattern.compile(pat,Pattern.CASE_INSENSITIVE);
		
		int[] no4sUIndices = EI.along(cond).filterInt(i->no4sUPattern.matcher(cond[i]).find()).toIntArray();
		boolean[] no4sU = new boolean[cond.length];
		for (int i : no4sUIndices)
			no4sU[i] = true;
		
		double[] polda = new double[cond.length];  Arrays.fill(polda, -1);
		double[] pnewa = new double[cond.length];  Arrays.fill(pnewa, -1);
		double[] dpolda = new double[cond.length];  Arrays.fill(dpolda, -1);
		double[] dpnewa = new double[cond.length];  Arrays.fill(dpnewa, -1);
		
		rit.sideEffect(a->a[0].startsWith("single_old"),a->parseRates(polda,a))
			.sideEffect(a->a[0].startsWith("single_new"),a->parseRates(pnewa,a))
			.sideEffect(a->a[0].startsWith("double_old"),a->parseRates(dpolda,a))
			.sideEffect(a->a[0].startsWith("double_new"),a->parseRates(dpnewa,a))
			.drain();
		
		HashMap<String,MutableInteger> g2l = new HashMap<>();
		genomic.getTranscripts().ei().forEachRemaining(tr->g2l.computeIfAbsent(tr.getData().getGeneId(), x->new MutableInteger()).max(tr.getRegion().getTotalLength())); 
		
		context.getLog().info("Preparing inference...");
		OptimNumericalIntegrationProportion[] vb = new OptimNumericalIntegrationProportion[cond.length];
		Progress pr = context.getProgress().init().setCount(vb.length);
		for (int i=0; i<vb.length; i++) {
			double pold = correct(polda[i]);
			double pnew = correct(pnewa[i]);
			double dpold = correct(dpolda[i]);
			double dpnew = correct(dpnewa[i]);
			if (doubleonly) {
				pold=pnew=0;
			}
			
			//context.getLog().info("For "+cond[i]+", using p_old = "+pold+"  p_new="+pnew+" dp_old="+dpold+" dp_new="+dpnew);
			
			vb[i] = no4sU[i]?null:new OptimNumericalIntegrationProportion(1, 1, pold,pnew,dpold,dpnew, lower, upper,true);
			int ui = i;
			pr.setDescription(()->cond[ui]).incrementProgress();
		}
		pr.finish();
		
		if (doubleonly) 
			context.getLog().info("Skipping non-overlapping parts of read pairs!");
		
		Function<String, String> sym = s->{
			if (s.endsWith("_intronic"))
				return genomic.getGeneTable("symbol").apply(StringUtils.removeFooter(s, "_intronic"))+"_intronic";
			return genomic.getGeneTable("symbol").apply(s);
		};
		
		context.getLog().info("Infering old/new proportion...");
		if (o10x)
			inferSparse(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, vb, data, getOutputFile(1));
		else
			infer(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, vb, data, getOutputFile(1));
		infer(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, vb, ext, getOutputFile(0));
		
		if (bias>0) {
			context.getLog().info("Infering old/new proportion for biased estimators...");
			double cbias = bias<1?1/bias:bias;
			OptimNumericalIntegrationProportion[] low = EI.wrap(vb).map(o->o.biasedEstimator(1, 1/cbias)).toArray(OptimNumericalIntegrationProportion.class);
			OptimNumericalIntegrationProportion[] high = EI.wrap(vb).map(o->o.biasedEstimator(1, cbias)).toArray(OptimNumericalIntegrationProportion.class);
			if (o10x) {
				inferSparse(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, low, data, new File(getOutputFile(1).getPath()+".low"));
				inferSparse(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, high, data, new File(getOutputFile(1).getPath()+".high"));
			}
			else {
				infer(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, low, data, new File(getOutputFile(1).getPath().replace(".tsv.gz", ".low.tsv.gz")));
				infer(context, sym, g2l, cond,no4sU, full, nthreads, lower, upper, high, data, new File(getOutputFile(1).getPath().replace(".tsv.gz", ".high.tsv.gz")));
			}
		}

		
		try {
			context.getLog().info("Running R scripts for plotting");
			RRunner r = new RRunner(prefix+".plotexonintron.R");
			r.set("prefix",prefix);
			r.addSource(getClass().getResourceAsStream("/resources/R/plotexin.R"));
			r.run(true);
		} catch (Throwable e) {
			context.getLog().log(Level.SEVERE, "Could not plot!", e);
		}
		
		return null;
	}


	private void inferSparse(GediProgramContext context, Function<String, String> sym, HashMap<String, MutableInteger> g2l, String[] cond, boolean[] no4sU, boolean full, int nthreads, double lower, double upper, OptimNumericalIntegrationProportion[] vb, File data, File dir) throws IOException {

		dir.mkdirs();
		if (!dir.isDirectory()) throw new IOException("Could not create directory "+dir);
		
		ArrayList<String> genes = new ArrayList<>();
		int nonzeros = 0;
		for (GeneData gd : EI.deserialize(new BinarySerializableSerializer<>(GeneData.class), new PageFile(data.getPath()), r->r.close()).loop()) {
			for (int c=0; c<gd.getReadCount().length(); c++)
				if (gd.getReadCount().getDouble(c)>0)
					nonzeros++;
			genes.add(gd.getGene()); 
		}
		
		int total = genes.size();
		
		EI.wrap(cond).print(new File(dir,"barcodes.tsv.gz").getPath());
		LineWriter out = new LineOrientedFile(dir, "features.tsv.gz").write();
		for (String g : genes) 
			out.writef("%s\t%s\tGene Expression\n", g,sym.apply(g));
		for (String g : genes) 
			out.writef("%s.NTR\t%s\tNTR\n", g,sym.apply(g));
		for (String g : genes) 
			out.writef("%s.CI\t%s\tCI\n", g,sym.apply(g));
		out.close();
		
		LineWriter writer = new LineOrientedFile(dir,"matrix.mtx.gz").write();
		writer.writeLine("%%MatrixMarket matrix coordinate real general");
		writer.writeLine("%metadata_json: {\"format_version\": 2, \"software_version\": \"3.0.2\"}");
		writer.writef("%d %d %d\n", genes.size()*3,cond.length,nonzeros*3);
		
		HashMap<String, Integer> geneIndex = EI.wrap(genes).indexPosition();
		
		EI.deserialize(new BinarySerializableSerializer<>(GeneData.class), new PageFile(data.getPath()), r->r.close())
			.progress(context.getProgress(), total, r->"Processing gene "+r.getGene())
			.parallelized(nthreads, 32, ei->ei.map(g->{
				return g.infer(vb,context.getLog());
			})
			.map(gp->{
				StringBuilder sb = new StringBuilder();
				for (int i=0; i<cond.length; i++) {
					if (gp.getReadCount().getDouble(i)>0) {
						SlamEstimationResult est = gp.getEstimated(i);
						int g = geneIndex.get(gp.getGene());
						sb.append(String.format("%d %d %d\n%d %d %.4f\n%d %d %.4f\n", 
									g+1,i+1,gp.getReadCount().getInt(i),
									g+genes.size()+1,i+1,est.getMap(),
									g+2*genes.size()+1,i+1,est.getUpper()-est.getLower()
								));
					}
				}
				sb.deleteCharAt(sb.length()-1);
				return sb.toString();
			})).print(writer);
		
		writer.close();
		

	}

	private void infer(GediProgramContext context, Function<String, String> sym, HashMap<String, MutableInteger> g2l, String[] cond, boolean[] no4sU, boolean full, int nthreads, double lower, double upper, OptimNumericalIntegrationProportion[] vb, File data, File outputFile) throws IOException {

		LineWriter writer = new LineOrientedFile(outputFile.getPath()).write();
		writer.write("Gene\tSymbol");
		for (int i=0; i<cond.length; i++) 
			writer.writef("\t%s Readcount",cond[i]);
		for (int i=0; i<cond.length; i++) 
			if (!no4sU[i])
				writer.writef("\t%s %.2f quantile\t%s Mean\t%s MAP\t%s %.2f quantile",cond[i],lower,cond[i],cond[i],cond[i],upper);
		for (int i=0; i<cond.length; i++) 
			if (!no4sU[i])
				writer.writef("\t%s alpha\t%s beta",cond[i],cond[i]);
		if (full)
			for (int i=0; i<cond.length; i++) 
				writer.writef("\t%s Conversions\t%s Coverage\t%s Double-Hits\t%s Double-Hit Coverage",cond[i],cond[i],cond[i],cond[i]);
		if (full)
			for (int i=0; i<cond.length; i++) 
				writer.writef("\t%s min2",cond[i]);
		writer.writeLine("\tLength");
		
		int total = (int) EI.deserialize(new BinarySerializableSerializer<>(GeneData.class), new PageFile(data.getPath()), r->r.close()).count();
		
		EI.deserialize(new BinarySerializableSerializer<>(GeneData.class), new PageFile(data.getPath()), r->r.close())
			.progress(context.getProgress(), total, r->"Processing gene "+r.getGene())
			.parallelized(nthreads, 32, ei->ei.map(g->{
				return g.infer(vb,context.getLog());
				
//				g.getConversionVector(0);
////				if (!g.getGene().equals("ERCC-00096"))
////					return null;
//				
//				SlamEstimationResult[] est = new SlamEstimationResult[cond.length];
//				
//				for (int i=0; i<cond.length; i++) {
//					est[i] = vb[i].infer(g.getConversionVector(i), g.getTotalVector(i), g.getWeightVector(i), g.getConversionDoubleVector(i), g.getTotalDoubleVector(i), g.getWeightDoubleVector(i), context.getLog());
//				}
//				
//				return new GeneProportion(g,est);
				
			})
			.map(gp->{
				StringBuilder sb = new StringBuilder();
				sb.append(gp.getGene()+"\t"+sym.apply(gp.getGene()));
				for (int i=0; i<cond.length; i++) 
					sb.append(String.format("\t%.04f",gp.getReadCount().getDouble(i)));
				for (int i=0; i<cond.length; i++) {
					if (!no4sU[i]) {
						SlamEstimationResult est = gp.getEstimated(i);
						sb.append(String.format("\t%.4f\t%.4f\t%.4f\t%.4f",est.getLower(),est.getMean(), est.getMap(),est.getUpper()));
					}
				}
				for (int i=0; i<cond.length; i++) { 
					if (!no4sU[i]) {
						double a = gp.getEstimated(i).getAlpha();
						double b = gp.getEstimated(i).getBeta();
						sb.append(String.format("\t%.4f\t%.4f",a,b));
					}
				}
				if (full) {
					for (int i=0; i<cond.length; i++) { 
						sb.append(String.format("\t%.4g\t%.4g\t%.4g\t%.4g",
								gp.getConversions().getDouble(i),
								gp.getCoverage().getDouble(i),
								gp.getDoubleHits().getDouble(i),
								gp.getDoubleHitCoverage().getDouble(i)
								));
					}
					NumericArray min2 = gp.getData().getReadsWithConversions(2);
					for (int i=0; i<cond.length; i++) { 
						sb.append(String.format("\t%.1f",min2.getDouble(i)));
					}
				}
				MutableInteger len = g2l.get(gp.getGene());
				sb.append("\t").append(len==null?0:len.N);
				return sb.toString();
			})).print(writer);
		writer.close();
	}



	private double correct(double d) {
		return Math.min(Math.max(d, 1E-6), 0.9);
	}



	private void parseRates(double[] rates, String[] a) {
		if (rates.length!=a.length-1) throw new RuntimeException("Rates file incompatible with reads!");
		for (int i=0; i<rates.length; i++)
			rates[i] = Double.parseDouble(a[i+1]);
	}

	
}
