package gedi.slam.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.regex.Pattern;

import javax.script.ScriptException;

import org.apache.commons.math3.stat.descriptive.rank.Median;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.slam.SlamParameterEstimation;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.tsv.formats.CsvReaderFactory;
import gedi.util.nashorn.JS;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;




/**
 * We need:
 * - the average sequencing error rate: 
 * 	Single-end (->Only one kind of reads): set to zero!
 * 	Paired-end: Mismatch rate - Doublehit rate (e.g. T->C mm in First read antisense - A->G double hit); robust fit over A->G and T->C
 * 
 * - the conversion (/transcription error/promiscuous editing) rate for old RNA:
 * 		- linear model trained from no4sU+IAA reads (with features T->A and T->G); take data from overlap (paired-end) or all
 * 
 * - the conversion (...) rate for new RNA:
 * 	Paired-end:
 * 		expected double hits from old = N*(conv_o + err^2); compute fraction among observed double hit read pairs f; this is a binom mix model with mixture weight and p1 known -> subtract expected from component 1 and estimate directly
 * 	Single end:
 * 		em and hope for the best! 
 * 
 * @author erhard
 *
 */
public class SlamRates extends GediProgram {

	
	public SlamRates(SlamParameterSet params) {
		addInput(params.mismatchFile);
		addInput(params.doublehitFile);
		addInput(params.reads);
		addInput(params.conv);
		addInput(params.err);
		addInput(params.no4sUpattern);
		addInput(params.binomFile);
		addInput(params.plot);
		addInput(params.prefix);
		addInput(params.minEstimateReads);
		addInput(params.errlm);
		addInput(params.errlm2);
		addInput(params.strandnessFile);
		addInput(params.dataFile);
		addInput(params.binomOverlapFile);
		
		addOutput(params.rateTable);
		addOutput(params.ntrFile);

	}
	
	
	public String execute(GediProgramContext context) throws IOException {
		
		File mmFile = getParameter(0);
		File dhFile = getParameter(1);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(2);
		String conv = getParameter(3);
		String err = getParameter(4);
		String pat = getParameter(5);
		File binom = getParameter(6);
		boolean plot = getBooleanParameter(7);
		String prefix = getParameter(8);
		int minEstimateReads = getIntParameter(9);
		String errlm = getParameter(10);
		String errlm2 = getParameter(11);
		Strandness strandness = Strandness.valueOf(EI.lines(this.<File>getParameter(12)).first());
		File data = getParameter(13);
		File binomOverlap = getParameter(14);
		boolean newMethod = conv!=null && conv.equals("new");
		
		if (strandness.equals(Strandness.AutoDetect))
			throw new RuntimeException("If you rerun, either delete the snps file or specify -strandness!");
		
		Pattern no4sUPattern = Pattern.compile(pat,Pattern.CASE_INSENSITIVE);
		
		String[] cond = reads.getMetaDataConditions();
		HashMap<String, Integer> cind = ArrayUtils.createIndexMap(cond);
		
		int[] no4sUIndices = EI.along(cond).filterInt(i->no4sUPattern.matcher(cond[i]).find()).toIntArray();
		boolean[] no4sU = new boolean[cond.length];
		for (int i : no4sUIndices)
			no4sU[i] = true;
		
		RateDataset rates = new RateDataset(cind, mmFile.getPath(), dhFile.getPath());
		
		if (no4sUIndices.length>=1) {
			context.getLog().info("Determining linear models using "+EI.wrap(no4sUIndices).map(i->cond[i]).concat(","));

			char[] errmod = {errlm.length()==2?errlm.charAt(0):'T',
							errlm.length()==2?errlm.charAt(1):'A'};
			
			
			double[] firstSense = EI.wrap(no4sUIndices).mapToDouble(i->Math.log(rates.getMismatchRate("Exonic", i, 'T', 'C',true,true)/rates.getMismatchRate("Exonic", i, errmod[0], errmod[1],true,true))).toDoubleArray();
			double[] firstAntisense = EI.wrap(no4sUIndices).mapToDouble(i->Math.log(rates.getMismatchRate("Exonic", i, 'A', 'G',true,false)/rates.getMismatchRate("Exonic", i, SequenceUtils.getDnaComplement(errmod[0]), SequenceUtils.getDnaComplement(errmod[1]),true,false))).toDoubleArray();
			double[] secondSense = EI.wrap(no4sUIndices).mapToDouble(i->Math.log(rates.getMismatchRate("Exonic", i, 'A', 'G',false,true)/rates.getMismatchRate("Exonic", i, SequenceUtils.getDnaComplement(errmod[0]), SequenceUtils.getDnaComplement(errmod[1]),false,true))).toDoubleArray();
			double[] secondAntisense = EI.wrap(no4sUIndices).mapToDouble(i->Math.log(rates.getMismatchRate("Exonic", i, 'T', 'C',false,false)/rates.getMismatchRate("Exonic", i, errmod[0], errmod[1],false,false))).toDoubleArray();
			double[] all = ArrayUtils.concat(firstSense,firstAntisense,secondSense,secondAntisense);
			
			
			
			double f = Math.exp(new Median().evaluate(all));
			
			if (Double.isNaN(f)) {
				context.getLog().warning("Could not learn linear model, no data! Falling back to default model.");
			}
			else {
				context.getLog().info("errlm="+String.valueOf(errmod)+"*"+f);
				errlm = String.valueOf(errmod)+"*"+f;
	
				if (!rates.isSingleEnd()) {
					f = Math.exp(new Median().evaluate(EI.wrap(no4sUIndices).mapToDouble(i->Math.log(rates.getDoubleHitRate("Exonic", i, 'T', 'C')/rates.getDoubleHitRate("Exonic", i, errmod[0], errmod[1]))).toDoubleArray()));
					context.getLog().info("errlm2="+String.valueOf(errmod)+"*"+f);
					errlm2 = String.valueOf(errmod)+"*"+f;
				}
			}
			
		}
		
		
		
		double[] single_old = new double[cond.length];
		
		if (StringUtils.isNumeric(err)) {
			Arrays.fill(single_old, Double.parseDouble(err));
		} else {
			
			
			if (strandness.equals(Strandness.Sense)) { 
				context.getLog().info("Estimate conversions for old RNA using linear model "+errlm+" for mismatches on exons");
				estimateOld(single_old,errlm, (c,g,r)->rates.getMismatchRate("Exonic",c,g,r,true,true));
			}
			else if (strandness.equals(Strandness.Antisense)) {
				context.getLog().info("Estimate conversions for old RNA using linear model "+errlm+" for mismatches on exons");
				estimateOld(single_old,errlm, (c,g,r)->rates.getMismatchRate("Exonic",c,SequenceUtils.getDnaComplement(g),SequenceUtils.getDnaComplement(r),true,false));
			} else {// unspecific 
				context.getLog().info("Estimate conversions for old RNA using linear model "+errlm+" for mismatches on exons");
				double[] a = new double[cond.length];
				estimateOld(single_old,errlm, (c,g,r)->rates.getMismatchRate("Exonic",c,g,r,true,true));
				estimateOld(a,errlm, (c,g,r)->rates.getMismatchRate("Exonic",c,SequenceUtils.getDnaComplement(g),SequenceUtils.getDnaComplement(r),true,false));
				for (int i=0; i<single_old.length; i++)
					single_old[i] = (single_old[i]+a[i])/2;
			}
		}
		
		double[] double_old = new double[cond.length];
		
		if (StringUtils.isNumeric(err)) {
			Arrays.fill(double_old, Double.parseDouble(err));
		} else {
			if (!rates.isSingleEnd()) {
				context.getLog().info("Estimate overlap conversions for old RNA using linear model "+errlm2+" for double hits on exons");
				estimateOld(double_old,errlm2, (c,g,r)->rates.getDoubleHitRate("Exonic",c,g,r));
			} 
		}
		
		
		
		
		
//		if (!rates.isSingleEnd()) {
//			context.getLog().info("Estimate conversions for new RNA overlap trick");
//			conv_new = SlamParameterEstimation.fromOverlapFile(overlap.getAbsolutePath(),r_err,conv_old,cind);
//		} else {
//			context.getLog().info("Estimate conversions for new RNA using EM");
//			SlamParameterEstimation est = new SlamParameterEstimation()
//					.setMinEstimateReads(minEstimateReads)
//					.setLogger(context.getLog());
//			conv_new = est.fromBinomFile(binom.getAbsolutePath(),getOutputFile(2).getPath(),r_err,conv_old);
//			didbinom = true;
//		}
		
		
		double[] single_new;
		double[] double_new;
		
		
		if (StringUtils.isNumeric(conv)) {
			single_new = new double[cond.length];
			double_new = new double[cond.length];
			Arrays.fill(single_new, Double.parseDouble(conv));
			Arrays.fill(double_new, Double.parseDouble(conv));
		} else {
			
		
			context.getLog().info("Estimate conversions for new RNA using "+(newMethod?"Nelder-Mead":"EM"));
			SlamParameterEstimation est = new SlamParameterEstimation(newMethod)
					.setMinEstimateReads(minEstimateReads)
					.setLogger(context.getLog());
			single_new = est.fromBinomFile(binom.getAbsolutePath(),cond,context.<SlamParameterSet>getParams().binomOut.getFile().getPath(),single_old,no4sU);
			
			double_new = new double[cond.length];
			if (!rates.isSingleEnd()) {
				context.getLog().info("Estimate conversions for new RNA in overlap using "+(newMethod?"Nelder-Mead":"EM"));
				est = new SlamParameterEstimation(newMethod)
						.setMinEstimateReads(minEstimateReads)
						.setLogger(context.getLog());
				double_new = est.fromBinomFile(binomOverlap.getAbsolutePath(),cond,context.<SlamParameterSet>getParams().binomOverlapOut.getFile().getPath(),double_old,no4sU);
			}
	
			if (rates.isSingleEnd() || double_new.length!=cond.length) 
				for (int i=0; i<cond.length; i++) 
					context.getLog().info(String.format("For %s, conv_old=%.3g conv_new=%.3g",cond[i],single_old[i],single_new[i]));
			else
				for (int i=0; i<cond.length; i++) 
					context.getLog().info(String.format("For %s, conv_old=%.3g conv_new=%.3g dconv_old=%.3g dconv_new=%.3g",cond[i],single_old[i],single_new[i],double_old[i],double_new[i]));
		
			LineWriter ntrOut = new LineOrientedFile(getOutputFile(1).getPath()).write().writeLine("Mode\tType\tCondition\tT->C\tp_new\tp_old\tntr_lower\tntr\tntr_upper");
			SlamParameterEstimation.writeNtrs(binom.getAbsolutePath(), "single",cond, ntrOut, single_old, single_new);
			SlamParameterEstimation.writeNtrs(binomOverlap.getAbsolutePath(), "double",cond, ntrOut, double_old, double_new);
			ntrOut.close();
		}
		
		LineWriter writer = new LineOrientedFile(getOutputFile(0).getPath()).write();
		writer.write("Rate");
		for (String c : reads.getMetaDataConditions())
			writer.write("\t"+c);
		writer.writeLine();
		writer.writef("single_old\t%s\n", NumericArray.wrap(single_old).toArrayString("\t",false));
		writer.writef("single_new\t%s\n", NumericArray.wrap(single_new).toArrayString("\t",false));
		if (!rates.isSingleEnd() && double_new.length==cond.length) {
			writer.writef("double_old\t%s\n", NumericArray.wrap(double_old).toArrayString("\t",false));
			writer.writef("double_new\t%s\n", NumericArray.wrap(double_new).toArrayString("\t",false));
		}
		writer.close();
		
//		HashMap<String, double[][]> tc = coll.get("TC");
		
//		Trie<double[]> named = new Trie<>();
		
		
//		// write out all T-C error rates for names
//		for (String name : tc.keySet()) {
//			writer.writef("%s\t%s\n", name,NumericArray.wrap(tc.get(name)).toArrayString("\t",false));
//			named.put(name.toLowerCase(),tc.get(name));
//		}
//		
//		// if there are no4sU samples: compute median T->C from Exonic
//		if (no4sUIndices.length>0) {
//			double[] rates = ArrayUtils.restrict(tc.get("Exonic"),no4sUIndices);
//			double median = new Median().evaluate(rates);
//			writer.writef("median_no4sU\t%s\n", EI.repeat(cond.length, ()->median).concat("\t"));
//			named.put("median_no4su",EI.repeat(cond.length, ()->median).toDoubleArray());
//		}
//				
//		// if lm given use Exonic to predict T->C error rate
//		if (errlm!=null) {
//			if (strandness.equals(Strandness.Unspecific))
//				throw new RuntimeException("Cannot evaluate linear model with unstranded data!");
//			writer.write("lm");
//			double[] lm = new double[cond.length];
//			for (int c=0; c<cond.length; c++) {
//				JS js = new JS();
//				DoubleArrayList rates = new DoubleArrayList();
//				for (String mm : coll.keySet()) {
//					double rate = coll.get(mm).get("Exonic")[c];
//					js.putVariable(mm, rate);
//					rates.add(rate);
//				}
//				double med = rates.evaluate(new Median());
//				js.putVariable("median", med);
//				try {
//					lm[c] = js.eval(errlm);
//					writer.writef("\t%s",lm[c]+"");
//				} catch (ScriptException e) {
//					throw new RuntimeException("Could not evalutate given lm: "+errlm,e);
//				}
//			}
//			named.put("lm", lm);
//			writer.writeLine();
//		}
//		
//		// if unstranded: use Exonic TC or AG from other strand
//		if (strandness.equals(Strandness.Unspecific)) {
//			writer.writef("antisense\t%s\n", NumericArray.wrap(coll.get("AG").get("Exonic")).toArrayString("\t",false));
//			named.put("antisense",coll.get("AG").get("Exonic"));
//		}
//		
//		// prioritize error rate approaches
//		// 1. Specified via -err
//		// 2. Specified via -errlm
//		// 3. ercc
//		// 4. antisense
//		// 5. median from no4sU
//		// 6. throw Exception
//		
//		double[] err_use = null;
//		if (err!=null) {
//			err_use = parseSpecifiedRate(err,named,cond.length);
//			context.getLog().info("Using specified error rates: "+StringUtils.toString(err_use));
//		}
//		else if (named.containsKey("lm")) {
//			err_use = named.get("lm");
//			context.getLog().info("Using linear model error rates: "+StringUtils.toString(err_use));
//		}
//		else if (named.containsKey("ercc")) {
//			err_use = named.get("ercc");
//			context.getLog().info("Using ERCC rates: "+StringUtils.toString(err_use));
//		}
//		else if (named.containsKey("antisense")) {
//			err_use = named.get("antisense");
//			context.getLog().info("Using antisense rates: "+StringUtils.toString(err_use));
//		}
//		else if (named.containsKey("median_no4su")) {
//			err_use = named.get("median_no4su");
//			context.getLog().info("Using median no4sU error rates: "+StringUtils.toString(err_use));
//		}
//		else {
//			writer.close();
//			throw new RuntimeException("Could not determine error rate. Either specify -err <1st col from rates file> or -errlm!");
//		}
//		
//		if (err_use.length!=cond.length) {
//			writer.close();
//			throw new RuntimeException("Invalid parameter -err!");
//		}
//			
//		writer.writef("err\t%s\n", NumericArray.wrap(err_use).toArrayString("\t",false));
//		
//		
//		context.getLog().info("Estimate conversions using EM");
//		SlamParameterEstimation est = new SlamParameterEstimation()
//				.setMinEstimateReads(minEstimateReads)
//				.setLogger(context.getLog());
//		double[] conv_em = est.fromBinomFile(binom.getAbsolutePath(),getOutputFile(1).getPath(),err_use);
//		writer.writef("conv_em\t%s\n", NumericArray.wrap(conv_em).toArrayString("\t",false));
//		
//		context.getLog().info("Estimate conversions using overlap");
//		double[] conv_overlap = SlamParameterEstimation.fromOverlapFile(overlap.getAbsolutePath());
//		double[] conv_overlap_n = SlamParameterEstimation.countsFromOverlapFile(overlap.getAbsolutePath());
//		if (conv_overlap.length>0) {
//			writer.writef("conv_overlap\t%s\n", NumericArray.wrap(conv_overlap).toArrayString("\t",false));
//			for (int i=0; i<conv_overlap.length; i++)
//				context.getLog().info(String.format("%s\tp_c=%.4g +/- %.4g", cond[i],conv_overlap[i], Math.sqrt(conv_overlap[i]*(1-conv_overlap[i])/conv_overlap_n[i] )  ));
//				
//		}
//		
//		
//		// prioritize conv rates
//		// 1. Specified via -conv
//		// 2. fit using overlap double hits
//		// 3. fit by EM on exonic with given error rates
//		double[] conv_use = null;
//		if (conv!=null) {
//			conv_use = parseSpecifiedRate(conv,named,cond.length);
//			context.getLog().info("Using specified conversion rates: "+StringUtils.toString(conv_use));
//		} else if (conv_overlap.length>0){
//			conv_use = conv_overlap;
//			context.getLog().info("Using conversion rates from overlap: "+StringUtils.toString(conv_use));
//		} else if (!est.isTooFew()){
//			conv_use = conv_em;
//			context.getLog().info("Using conversion rates from EM algorithm: "+StringUtils.toString(conv_use));
//		} else {
//			writer.close();
//			throw new RuntimeException("Could not determine conversion rate. Specify -conv <1st col from rates file>!");
//		}
//		
//		writer.writef("conv\t%s\n", NumericArray.wrap(conv_use).toArrayString("\t",false));
//		
//		
//		writer.close();
//		
//		
		if (plot) {
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".plotbinom.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/plotbinom.R"));
				r.run(true);

				context.getLog().info("Running R scripts for plotting");
				r = new RRunner(prefix+".plotrates.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/plotrates.R"));
				r.run(true);
				
//				context.getLog().info("Running R scripts for plotting");
//				r = new RRunner(prefix+".ploterrorate.R");
//				r.set("prefix",prefix);
//				r.addSource(getClass().getResourceAsStream("/resources/R/ploterrorate.R"));
//				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		
		return null;
	}


	private double estimateError(RateDataset rates, String path, String[] cond) throws IOException {
		// Mismatch rate - Doublehit rate (e.g. T->C mm in First read antisense - A->G double hit); robust fit over A->G and T->C
		
		DoubleArrayList[] mm = new DoubleArrayList[4];
		DoubleArrayList[] dh = new DoubleArrayList[4];
		for (int i=0; i<4; i++) {
			mm[i] = new DoubleArrayList();
			dh[i] = new DoubleArrayList();
		}
		
		LineWriter writer = new LineOrientedFile(path).write();
		writer.writeLine("Condition\tType\tMismatch\tDoublehit");
		for (int c=0; c<rates.getNumConditions(); c++) {
			mm[0].add( (rates.getMismatches("Exonic", c, 'T', 'C', true, false)+rates.getMismatches("Exonic", c, 'T', 'C', false, true))
					/ (rates.getMismatchCoverage("Exonic", c, 'T', 'C', true, false)+rates.getMismatchCoverage("Exonic", c, 'T', 'C', false, true))
					);
			dh[0].add(rates.getDoubleHitRate("Exonic", c, 'A', 'G'));
			writer.writef("%s\t%s\t%.5g\t%.5g\n",cond[c],"antisense T->C",mm[0].getLastDouble(),dh[0].getLastDouble());
			
			mm[1].add( (rates.getMismatches("Exonic", c, 'A', 'G', true, true)+rates.getMismatches("Exonic", c, 'A', 'G', false, false))
					/ (rates.getMismatchCoverage("Exonic", c, 'A', 'G', true, true)+rates.getMismatchCoverage("Exonic", c, 'A', 'G', false, false))
					);
			dh[1].add(rates.getDoubleHitRate("Exonic", c, 'A', 'G'));
			writer.writef("%s\t%s\t%.5g\t%.5g\n",cond[c],"sense A->G",mm[1].getLastDouble(),dh[1].getLastDouble());
			
			
			mm[2].add( (rates.getMismatches("Exonic", c, 'T', 'C', true, true)+rates.getMismatches("Exonic", c, 'T', 'C', false, false))
					/ (rates.getMismatchCoverage("Exonic", c, 'T', 'C', true, true)+rates.getMismatchCoverage("Exonic", c, 'T', 'C', false, false))
					);
			dh[2].add(rates.getDoubleHitRate("Exonic", c, 'T', 'C'));
			writer.writef("%s\t%s\t%.5g\t%.5g\n",cond[c],"sense T->C",mm[2].getLastDouble(),dh[2].getLastDouble());
			
			mm[3].add( (rates.getMismatches("Exonic", c, 'A', 'G', true, false)+rates.getMismatches("Exonic", c, 'A', 'G', false, true))
					/ (rates.getMismatchCoverage("Exonic", c, 'A', 'G', true, false)+rates.getMismatchCoverage("Exonic", c, 'A', 'G', false, true))
					);
			dh[3].add(rates.getDoubleHitRate("Exonic", c, 'T', 'C'));
			writer.writef("%s\t%s\t%.5g\t%.5g\n",cond[c],"antisense A->G",mm[3].getLastDouble(),dh[3].getLastDouble());
		}
		writer.close();
		
		double sum = 0;
		for (int i=0; i<4; i++) {
			NumericArray co = mm[i].toNumericArray().copy();
			co.subtract(dh[i].toNumericArray());
			sum+=co.evaluate(NumericArrayFunction.Median);
		}
		
		return sum/4;
	}

	private void estimateOld(double[] re, String errlm, IntCharCharToDoubleFunction fun) {
		char[] nucl = {'A','C','G','T'};
		for (int c=0; c<re.length; c++) {
			JS js = new JS();
			DoubleArrayList rates = new DoubleArrayList();
			for (char g : nucl)for (char r : nucl) if (g!=r){
				double rate = fun.apply(c, g, r);
				js.putVariable(String.valueOf(new char[] {g,r}), rate);
				rates.add(rate);
			}
			double med = rates.evaluate(new Median());
			js.putVariable("median", med);
			try {
				re[c] = js.eval(errlm);
			} catch (ScriptException e) {
				throw new RuntimeException("Could not evalutate given lm: "+errlm,e);
			}
		}
	}

	private double[] parseSpecifiedRate(String err, Trie<double[]> named, int cond) {
		double[] re = named.getUniqueWithPrefix(err.toLowerCase());
		if (re!=null) return re;
		re = EI.split(err, ',').mapToDouble(s->Double.parseDouble(s)).toDoubleArray();
		if(re.length==1)
			re = EI.repeat(cond,re[0]).toDoubleArray();
		return re;
	}

	@FunctionalInterface
	private static interface IntCharCharToDoubleFunction {
		double apply(int i, char g, char r);
	}
	
	
	public static class RateDataset {
		
		
		HashMap<String,double[][][][][][]> mismatches = new HashMap<>();
		HashMap<String,double[][][][]> doubles = new HashMap<>();;
		
		public RateDataset(HashMap<String,Integer> cind, String mm, String dh) {
			for (MismatchLine l : new CsvReaderFactory().setParseToClass(MismatchLine.class).createReader(mm).iterateObjects(MismatchLine.class).loop()) {
				boolean sense = l.Category.endsWith("Sense");
				boolean antisense = l.Category.endsWith("Antisense");
				String cat = StringUtils.removeFooter(StringUtils.removeFooter(l.Category, "Sense"),"Antisense");
				int cond = cind.get(l.Condition);
				boolean first = l.Orientation.equals("First");
				int g = SequenceUtils.inv_nucleotides[l.Genomic.charAt(0)];
				int r = SequenceUtils.inv_nucleotides[l.Read.charAt(0)];
				if (sense||antisense) {
					double[][][][][][] a = mismatches.computeIfAbsent(cat, x->new double[4][4][2][2][2][cind.size()]);
					a[g][r][first?0:1][sense?0:1][0][cond] = l.Mismatches;
					a[g][r][first?0:1][sense?0:1][1][cond] = l.Coverage;
				}
				
			}
			
			for (DoubleHitLine l : new CsvReaderFactory().setParseToClass(DoubleHitLine.class).createReader(dh).iterateObjects(DoubleHitLine.class).loop()) {
				boolean sense = l.Category.endsWith("Sense");
				boolean antisense = l.Category.endsWith("Antisense");
				String cat = StringUtils.removeFooter(StringUtils.removeFooter(l.Category, "Sense"),"Antisense");
				int cond = cind.get(l.Condition);
				int g = SequenceUtils.inv_nucleotides[l.Genomic.charAt(0)];
				int r = SequenceUtils.inv_nucleotides[l.Read.charAt(0)];
				if (!sense&&!antisense) {
					double[][][][] a = doubles.computeIfAbsent(cat, x->new double[4][4][2][cind.size()]);
					a[g][r][0][cond] = l.Hits;
					a[g][r][1][cond] = l.Coverage;
				}
				
			}
			
		}
		
		public int getNumConditions() {
			return mismatches.get("Exonic")[0][1][0][0][0].length;
		}

		public boolean isSingleEnd() {
			return ArrayUtils.sum(mismatches.get("Exonic")[0][1][1][0][1])+ArrayUtils.sum(mismatches.get("Exonic")[0][1][1][1][1])==0;
		}
		
		public double getMismatchRate(String category, int cond, char genomic, char read, boolean first, boolean sense) {
			double[][] a = mismatches.get(category)[SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][first?0:1][sense?0:1];
			return a[0][cond]/a[1][cond];
		}
		
		public double getMismatchCoverage(String category, int cond, char genomic, char read, boolean first, boolean sense) {
			double[][] a = mismatches.get(category)[SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][first?0:1][sense?0:1];
			return a[1][cond];
		}
		
		public double getMismatches(String category, int cond, char genomic, char read, boolean first, boolean sense) {
			double[][] a = mismatches.get(category)[SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][first?0:1][sense?0:1];
			return a[0][cond];
		}
		
		public double getDoubleHitRate(String category, int cond, char genomic, char read) {
			double[][] a = doubles.get(category)[SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]];
			return a[0][cond]/a[1][cond];

		}
		
		
	}
	
	private static class MismatchLine {
		private String Category;
		private String Condition;
		private String Orientation;
		private String Genomic;
		private String Read;
		private long Coverage;
		private long Mismatches;
	}
	
	private static class DoubleHitLine {
		private String Category;
		private String Condition;
		private String Genomic;
		private String Read;
		private long Coverage;
		private long Hits;
	}
	
	
}

