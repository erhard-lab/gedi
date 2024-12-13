package gedi.slam;

import static java.lang.Math.exp;
import static java.lang.Math.log;
import static java.lang.Math.max;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.Predicate;
import java.util.logging.Logger;

import org.apache.commons.math3.exception.TooManyEvaluationsException;

import gedi.grand3.estimation.models.Grand3Model;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.KahanSummation;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import jdistlib.Binomial;
import jdistlib.math.MultivariableFunction;
import jdistlib.math.opt.Bobyqa;
import jdistlib.math.opt.BobyqaConfig;
import jdistlib.math.opt.OptimizationResult;

public class SlamParameterEstimation {

	
	private double maxErrComponentFraction;
	private double precision = 1E-12;
	private int minTotal;
	
	private Logger log;
	private boolean toofew = false;
	
	private boolean newMethod = false;
	
	public SlamParameterEstimation(boolean newMethod) {
		this(newMethod,10000,0.01);
	}
	
	public SlamParameterEstimation(boolean newMethod, int minTotal, double maxErrComponentFraction) {
		super();
		this.newMethod = newMethod;
		this.minTotal = minTotal;
		this.maxErrComponentFraction = maxErrComponentFraction;
	}
	
	public SlamParameterEstimation setMinEstimateReads(int reads) {
		this.minTotal = reads;
		return this;
	}
	
	public SlamParameterEstimation setLogger(Logger log) {
		this.log = log;
		return this;
	}
	
	public boolean isTooFew() {
		return toofew;
	}
	
	
	private static final double lse(double u, double v) {
	    double m = max(u,v);
	    return log(exp(u-m)+exp(v-m))+m;
	}
	public double estimateNew(double[][] x, double errorProb, boolean no4sU) {
		for (int i=0; i<x.length; i++)
			if (x[0].length!=x[i].length) throw new RuntimeException("Not a matrix!");
		
		double e = errorProb;

		MultivariableFunction fun = a->{
			double m = a[0];
			double c = a[1];
			
			double sum = 0;
			for (int n=0; n<x[0].length; n++) 
				for (int k=0; k<x.length; k++)
					if (x[k][n]>0)
						sum+=lse(Math.log(m)+Binomial.density(k, n, e, true), Math.log(1-m)+Binomial.density(k, n, c, true))*x[k][n];
				
			return sum;
		};
		BobyqaConfig bcfg = new BobyqaConfig(new double[] {0.5,0.04},new double[] {0,0},new double[] {1,1},fun,100_000,false);
		bcfg.setNumInterpolationPoints(5);
		
		OptimizationResult res = Bobyqa.bobyqa(bcfg.getInitialGuess(), bcfg.getLowerBound(), bcfg.getUpperBound(),
				bcfg.getObjectiveFunction(), bcfg.getNumInterpolationPoints(), bcfg.getInitialTrustRegionRadius(),
				bcfg.getStoppingTrustRegionRadius(), bcfg.getMaxNumFunctionCall(), bcfg.isMinimize());
		
		
        double[] mc = res.mX;
		
//        if (!nm.getConvStatus() && log!=null && !no4sU) {
//			log.warning("Did not converge; maybe your experiments did not work!");
//		}
        
        
////		try {
//		double[] mc = new SimplexOptimizer(1E-5,1E-5).optimize(
//				GoalType.MAXIMIZE,
//				new NelderMeadSimplex(2),
//				new ObjectiveFunction(fun),
//				new InitialGuess(new double[]{0.5,0.04}),
//				new MaxEval(10000)
//				).getPointRef();
		
		return mc[1];
		
		
	}
	
	
	/**
	 * re is {err,conv,p_l,p,p_h}
	 * @param x
	 * @return
	 * @throws EstimationException 
	 */
	public double estimate(double[][] x, double errorProb, boolean no4sU) {
		if (newMethod) {
			return estimateNew(x, errorProb,no4sU);
		}
		
		
		for (int i=0; i<x.length; i++)
			if (x[0].length!=x[i].length) throw new RuntimeException("Not a matrix!");
		
		
		boolean[][] oc = computeOnlyConversionComponent(x,errorProb);
		double total = getTotalReads(x,oc);
		toofew = total<minTotal;
		if (toofew && log!=null && !no4sU) {
			log.warning("Too few reads ("+(int)total+") usable for estimating error and conversion rates; maybe your experiments did not work. To circumvent this warning, invoke with smaller -minEstimateReads or specify rates directly via -conv/-err!");
		}
		
		double[][] work = ArrayUtils.cloneMatrix(x);
		
		double conv = estimateConversion(work,oc);
		
		
//		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
//		GrandSlamMaskedBinomialModel model = new GrandSlamMaskedBinomialModel();
//		MaximumLikelihoodParametrization re = estimator.fit(model,model.createPar(0.05), new MutablePair<>(x,oc));
//		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(model, re, new MutablePair<>(x,oc));
//		
//		System.out.println(re);
//
//		GrandSlamMaskedTBinomialModel model2 = new GrandSlamMaskedTBinomialModel();
//		MaximumLikelihoodParametrization re2 = estimator.fit(model2,model2.createPar(0.05,4), new MutableTriple<>(x,oc,errorProb));
//		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(model2, re2, new MutableTriple<>(x,oc,errorProb));
//		System.out.println(re2);

		return conv;
	}
	


	private static class GrandSlamMaskedBinomialModel implements MaximumLikelihoodModel<MutablePair<double[][],boolean[][]>> {
		
		
		public double logLik(MutablePair<double[][],boolean[][]> data, double[] param) {
			double[][] x = data.Item1;
			boolean[][] oc = data.Item2;
			double p = param[0];
			
			double re = 0;
			for (int n=0; n<oc[0].length; n++) {
				double norm = 0;
				double sum = 0;
				for (int k=0; k<=Math.min(n, oc.length-1); k++) 
					if (oc[k][n]) {
						re+=Grand3Model.lbinom(k, n, p)*x[k][n];
						sum+=x[k][n];
						norm=norm==0?Grand3Model.lbinom(k, n, p):Grand3Model.lse(norm, Grand3Model.lbinom(k, n, p));
					}
				if (sum>0) 
					re-=norm*sum;
			}
			return re;
		}

		
		public MaximumLikelihoodParametrization createPar() {
			return createPar(0.05);
		}
	
		public MaximumLikelihoodParametrization createPar(double pconv) {
			return new MaximumLikelihoodParametrization(
					"maskedbinom",
					new String[] {"p.conv"},
					new double[] {pconv},
					new String[] {"%.4f"},
					new double[] {0},
					new double[] {1}
				);
		}
	}


	private static class GrandSlamMaskedTBinomialModel implements MaximumLikelihoodModel<MutableTriple<double[][],boolean[][],Double>> {
		
		
		public double logLik(MutableTriple<double[][],boolean[][],Double> data, double[] param) {
			double[][] x = data.Item1;
			boolean[][] oc = data.Item2;
			double perr=data.Item3;
			
			double p = param[0];
			double shape = param[1];
			
			double re = 0;
			for (int n=0; n<oc[0].length; n++) {
				double norm = 0;
				double sum = 0;
				for (int k=0; k<=Math.min(n, oc.length-1); k++) 
					if (oc[k][n]) {
						re+=Grand3Model.ldtbbinom(k, n, perr,p,shape)*x[k][n];
						sum+=x[k][n];
						norm=norm==0?Grand3Model.ldtbbinom(k, n, perr,p,shape):Grand3Model.lse(norm, Grand3Model.ldtbbinom(k, n, perr,p,shape));
					}
				if (sum>0) 
					re-=norm*sum;
			}
			return re;
		}

		
		public MaximumLikelihoodParametrization createPar() {
			return createPar(0.05,5);
		}
	
		public MaximumLikelihoodParametrization createPar(double pconv, double shape) {
			return new MaximumLikelihoodParametrization(
					"maskedbinom",
					new String[] {"p.conv","shape"},
					new double[] {pconv,shape},
					new String[] {"%.4f","%.2f"},
					new double[] {0,-2},
					new double[] {1,5}
				);
		}
	}


	private double getTotalReads(double[][] x, boolean[][] oc) {
		double re = 0;
		for (int n=0; n<oc[0].length; n++) 
			for (int k=0; k<oc.length; k++)
				if (oc[k][n])
					re+=x[k][n];
		return re;
	}

	private double estimateConversion(double[][] x, boolean[][] oc) {
		double low = 0;
		double high = 1;
		while (high-low>precision) {
			double p = (low+high)/2;
			eStep(x,p,oc);
			double np = mStep(x);
//			System.out.println(String.format("%.5g\t%.5g\t%.5g\t%.5g",low,high,p,np));
			
			if (np<p)
				high = p;
			else
				low = p;
		}
		return (low+high)/2;
	}

	private boolean[][] computeOnlyConversionComponent(double[][] x,double errorProb) {
		boolean[][] oc = new boolean[x.length][x[0].length];
		for (int n=0; n<oc[0].length; n++) {
			double s=0;
			for (int k=0; k<oc.length; k++) 
				s+=x[k][n];
			for (int k=0; k<oc.length; k++) {
				double expect = Binomial.density(k, n, errorProb, false)*s;//new BinomialDistribution(n, maxErrorProb).probability(k)*s;//Binomial.density(k, n, maxErrorProb, false)*s;
				oc[k][n] = expect<maxErrComponentFraction*x[k][n];
				if (k>0)
					oc[k][n] |= oc[k-1][n];
			}
			// only columns with at least two non-zero entries are useful!
			int using = 0;
			for (int k=0; k<oc.length; k++) {
				if (x[k][n]>0 && oc[k][n])
					using++;
			}
			if (using<2) {
				for (int k=0; k<oc.length; k++)
					oc[k][n]=false;
			}
		}
		return oc;
	}
	
	private boolean[][] computeOnlyZeros(double[][] x) {
		boolean[][] oc = new boolean[x.length][x[0].length];
		for (int n=0; n<oc[0].length; n++) 
			for (int k=0; k<oc.length; k++) 
				oc[k][n] = k==0;
		return oc;
	}

	private void eStep(double[][] x, double p, boolean[][] oc) {
		
		for (int n=0; n<x[0].length; n++) {
			double sx = 0;
			double sp = 0;
			for (int k=0; k<x.length; k++) {
				if (oc[k][n]) {
					sx+=x[k][n];
					sp+=Binomial.density(k, n, p, false);//new BinomialDistribution(n, p).probability(k);//Binomial.density(k, n, p, false);
				}
			}
			
			for (int k=0; k<x.length; k++) 
				if (!oc[k][n]) {
					x[k][n]=sx/sp*Binomial.density(k, n, p, false);//new BinomialDistribution(n, p).probability(k);//Binomial.density(k, n, p, false);
					if (Double.isNaN(x[k][n]))
						x[k][n] = 0;
				}
		}
	}
	
	private boolean isPToMuch(double[][] x, double[][] init, double p, boolean[][] oc) {
		
		for (int n=0; n<x[0].length; n++) {
			double sx = 0;
			double sp = 0;
			for (int k=0; k<x.length; k++) {
				if (oc[k][n]) {
					sx+=x[k][n];
					sp+=Binomial.density(k, n, p, false);//new BinomialDistribution(n, p).probability(k);//Binomial.density(k, n, p, false);
				}
			}
			
			for (int k=0; k<x.length; k++) 
				if (!oc[k][n]) {
					if (init[k][n]<sx/sp*Binomial.density(k, n, p, false))//new BinomialDistribution(n, p).probability(k))
						return true;
				}
		}
		return false;
	}

	private static double mStep(double[][] x) {
		KahanSummation num = new KahanSummation();
		KahanSummation den = new KahanSummation();
		for (int k=0; k<x.length; k++) {
			KahanSummation s = new KahanSummation();
			for (int n=0; n<x[0].length; n++) {
				s.add(x[k][n]);
			}
			num.add(k*s.getSum());
		}
		
		for (int n=0; n<x[0].length; n++) {
			KahanSummation s = new KahanSummation();
			for (int k=0; k<x.length; k++) {
				s.add(x[k][n]);
			}
			den.add(n*s.getSum());
		}
//		System.out.println(num+" "+den);
		return num.getSum()/den.getSum();
	}
	
	private double mStep(double[][] x, double[][] m) {
		KahanSummation num = new KahanSummation();
		KahanSummation den = new KahanSummation();
		for (int k=0; k<x.length; k++) {
			KahanSummation s = new KahanSummation();
			for (int n=0; n<x[0].length; n++) {
				s.add(m[k][n]*x[k][n]);
			}
			num.add(k*s.getSum());
		}
		
		for (int n=0; n<x[0].length; n++) {
			KahanSummation s = new KahanSummation();
			for (int k=0; k<x.length; k++) {
				s.add(m[k][n]*x[k][n]);
			}
			den.add(n*s.getSum());
		}
//		System.out.println(num+" "+den);
		return num.getSum()/den.getSum();
	}
	
	public static BinomialMatrix[] readBinomFile(String path, String type, String[] conditions) throws IOException {
		String[][] t = EI.lines(path).map(s->StringUtils.split(s, '\t')).toArray(String[].class);
		
		int[] cols = {0,1,2,3};
		Predicate<String[]> check = a->true;
		
		if (t.length==0)
			throw new RuntimeException("Invalid file format!");
		if (Arrays.equals(t[0],new String[] {"n","d","Condition","Type","count"})) {
			cols = new int[] {0,1,2,4};
			check = a->a[3].equals(type);
		} 
		else if (!Arrays.equals(t[0],new String[] {"n","d","Condition","count"}))
			throw new RuntimeException("Invalid file format!");
		if (t.length==1) return new BinomialMatrix[0];
		
		t = ArrayUtils.slice(t, 1);
		
		HashMap<String, Integer> cind = ArrayUtils.createIndexMap(conditions);
		
		int k = EI.wrap(t).mapToInt(l->Integer.parseInt(l[1])).max()+1;
		int n = EI.wrap(t).mapToInt(l->Integer.parseInt(l[0])).max()+1;
		
		
		BinomialMatrix[] xx = new BinomialMatrix[conditions.length];
		for (int i=0; i<xx.length; i++)
			xx[i] = new BinomialMatrix(k,n, conditions[i]);
		for (int l=0; l<t.length; l++) 
			if (check.test(t[l])) {
				int hcond = cind.get(t[l][cols[2]]);
				int hn = Integer.parseInt(t[l][cols[0]]);
				int hk = Integer.parseInt(t[l][cols[1]]);
				double hx = Double.parseDouble(t[l][cols[3]]);
				xx[hcond].cname = t[l][cols[2]];
				xx[hcond].A[hk][hn] += hx;
			}
		
		return xx;
	}
	
	public static class BinomialMatrix {
		String cname;
		private double[][] A;
		public BinomialMatrix(int k, int n, String cname) {
			A = new double[k][n];
			this.cname = cname;
		}
		
		public double estimate() {
			return mStep(A);
		}
		
		public double conversions() {
			KahanSummation re = new KahanSummation();
			for (int n=0; n<A[0].length; n++) {
				for (int k=1; k<A.length; k++) {
					re.add(k*A[k][n]);
				}
			}
			return re.getSum();
		}
		
		public double total() {
			KahanSummation re = new KahanSummation();
			for (int n=0; n<A[0].length; n++) {
				KahanSummation s = new KahanSummation();
				for (int k=0; k<A.length; k++) {
					s.add(A[k][n]);
				}
				re.add(n*s.getSum());
			}
			return re.getSum();
		}

		public double estimate(double err, double conv_old) {
			//expected double hits from old = N*(conv_o + err^2); compute fraction among observed double hit read pairs f; 
			// this is a binom mix model with mixture weight and p1 known -> subtract expected from component 1 and estimate directly
			
			double[][] As = ArrayUtils.cloneMatrix(A);
			double exp = ArrayUtils.sum(A)*(conv_old+err*err);
			double f = exp/total();
			for (int n=0; n<A[0].length; n++) {
				int N = 0;
				for (int k=0; k<A.length; k++) 
					N+=A[k][n];
					
				for (int k=0; k<A.length; k++) 
					As[k][n]=Math.max(0, As[k][n]-N*f*Binomial.density(k, n, conv_old, false));
			}
			
			return mStep(As);
		}

		public double getMismatchFrequency() {
			return conversions()/total();
		}
		
	}
	
//	public static double[] fromOverlapFile(String path, double r_err, double[] conv_old, HashMap<String,Integer> cind) throws IOException {
//		return EI.wrap(readBinomFile(path)).mapToDouble(m->m.estimate(r_err,conv_old[cind.get(m.cname)])).toDoubleArray();
//	}
//	
//	public static double[] countsFromOverlapFile(String path) throws IOException {
//		return EI.wrap(readBinomFile(path)).mapToDouble(m->m.total()).toDoubleArray();
//	}
	
	public static void main(String[] args) throws IOException {
		SlamParameterEstimation est = new SlamParameterEstimation(true);
		BinomialMatrix[] As = readBinomFile("/home/erhard/tmp/lior/run_with_mereged/slam.binom.tsv", "", new String[] {"2_D_alk_75","4sUTP_0_hpf_75","4sUTP_0_hpf_comb_75","uninj_no_alk_75_forbg"});
		System.out.println(As.length);
		System.out.println(As[2].A.length);
		double pre = Math.max(0, est.estimate(As[2].A,0.0006578,false));
		System.out.println(Arrays.toString(As[2].A[0]));
		System.out.println(Arrays.toString(As[2].A[1]));
		System.out.println(Arrays.toString(As[2].A[51]));
		System.out.println(Arrays.toString(As[2].A[52]));
		System.out.println(pre);
	}
	
	
	/**
	 * re is [conditions]
	 * @param path
	 * @return
	 * @throws IOException
	 */
	public double[] fromBinomFile(String path, String[] conditions, String output, double[] conv_old, boolean[] no4sU) throws IOException {
		
		BinomialMatrix[] As = readBinomFile(path, "", conditions);
		int cond = As.length;
		
		LineWriter out = output!=null?new LineOrientedFile(output).write().writeLine("Condition\tconv_old\tconv_new\tp_new_lower\tp_new\tp_new_upper\tStatus"):null;

		
		double[] re = new double[cond];
		for (int i=0; i<re.length; i++) {
			re[i] = Math.max(0, estimate(As[i].A,conv_old[i],no4sU[i]));
			SlamEstimationResult p = new SlamEstimationResult(0, 0.5, 0.5, 1);
			try {
				p = new OptimNumericalIntegrationProportion(1, 1, conv_old[i], re[i], conv_old[i], re[i], 0.05, 0.95, true).infer(As[i].A, null);
			} catch (TooManyEvaluationsException e) {
				if (log!=null)
					log.warning("Could not compute new fraction!");
			}
			if (log!=null)
				log.info("For "+As[i].cname+": old="+conv_old[i]+" new="+re[i]+" p in ["+p.getLower()+","+p.getUpper()+"]");
			if (out!=null)
				out.writef("%s\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%s\n", As[i].cname, conv_old[i], re[i], p.getLower(), p.getMean(), p.getUpper(), isTooFew()?"Too few":"Ok");
			
//			try {
//				FileUtils.writeAllText(ArrayUtils.matrixToString(xx[i]),new File(cnames[i]+"_test.mat"));
//			} catch (IOException e) {
//			}
			
			if (toofew && re[i]>0.1)
				re[i] = 0.1;
		}
		
		if (out!=null)
			out.close();
		
		return re;
	}

	
	public static void writeNtrs(String path, String mode, String[] conditions, LineWriter out, double[] conv_old, double[] conv_new) throws IOException {
		
		HeaderLine h = new HeaderLine();
		String[] types = EI.lines(path).header(h).split('\t').map(a->a[h.get("Type")]).unique(false).filter(s->s.length()>0).toArray(String.class);
		
		for (String type : types) {
			BinomialMatrix[] As = readBinomFile(path, type, conditions);
			
			for (int i=0; i<As.length; i++) {
				SlamEstimationResult p = new SlamEstimationResult(0, 0.5, 0.5, 1);
				try {
					p = new OptimNumericalIntegrationProportion(1, 1, conv_old[i], conv_new[i], conv_old[i], conv_new[i], 0.05, 0.95, true).infer(As[i].A, null);
				} catch (TooManyEvaluationsException e) {
				}
				if (out!=null)
					out.writef("%s\t%s\t%s\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n", mode, type, As[i].cname, 
							As[i].getMismatchFrequency(),conv_new[i],conv_old[i],
							p.getLower(), p.getMap(), p.getUpper());
			}
		}
		
	}

	
	/**
	 * re is [conditions][{err,conv,p}]
	 * @param path
	 * @return
	 * @throws IOException
	 */
	public double[] fromGeneData(GeneData[] a, double[] errorProbs) {
		
		int cond = a[0].getReads()[0].getCount().length();
		int n = EI.wrap(a).unfold(gd->EI.wrap(gd.getReads())).mapToInt(rd->rd.getTotal()).max()+1;
		int k = EI.wrap(a).unfold(gd->EI.wrap(gd.getReads())).mapToInt(rd->rd.getConversions()).max()+1;
		
		double[][][] xx = new double[cond][k][n];
		
		for (ReadData d : EI.wrap(a).unfold(gd->EI.wrap(gd.getReads())).loop()) {
			for (int hcond=0; hcond<d.getCount().length(); hcond++) {
				int hn = d.getTotal();
				int hk = d.getConversions();
				double hx = d.getCount().get(hcond);
				xx[hcond][hk][hn] += hx;
			}
		}
		
//		try {
//			FileUtils.writeAllText(ArrayUtils.matrixToString(xx[0]),new File("test.mat"));
//		} catch (IOException e) {
//		}
		
		double[] re = new double[cond];
		for (int i=0; i<re.length; i++) {
			
			re[i] = estimate(xx[i],errorProbs[i],false);
		}
		
		return re;
		
	}
	
	
	public double[] fromText(String path, double[] errorProbs) throws IOException {
		double[][] m = EI.lines(path).map(s->StringUtils.parseDouble(StringUtils.split(s, '\t'))).toArray(double[].class);
		double[][][] xx = {m};
		double[] re = new double[1];
		for (int i=0; i<re.length; i++) {
			
			re[i] = estimate(xx[i],errorProbs[i],false);
		}
		
		return re;
		
	}
	
	
}
