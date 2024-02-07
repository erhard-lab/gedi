package gedi.grand3.estimation.models;

import static jdistlib.math.MathFunctions.lchoose;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gedi.grand3.estimation.IncompleteBeta;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.grand3.knmatrix.KnMatrixKey;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.functions.EI;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import jdistlib.Beta;
import jdistlib.Binomial;
import jdistlib.Binomial.RandomState;
import jdistlib.Uniform;
import jdistlib.math.MathFunctions;
import jdistlib.rng.RandomEngine;
import net.jafama.FastMath;

public abstract class Grand3Model implements MaximumLikelihoodModel<KNMatrixElement[]> {
	
	public static final double logdiff(double u, double v) {
		return u+FastMath.log1p(-FastMath.exp(-(u-v)));
	}
	
	/**
	 * The elki version is much faster!
	 * @param u
	 * @param v
	 * @return
	 */
	public static final double lse(double u, double v) {
		double max = Math.max(u, v);
		final double cutoff = max - 35.350506209; // log_e(2**51)
	    double acc = 0.;
	      if(v > cutoff) {
	        acc += v < max ? FastMath.exp(v - max) : 1.;
	      }
	      if(u > cutoff) {
		        acc += u < max ? FastMath.exp(u - max) : 1.;
		      }
	    return acc > 1. ? (max + FastMath.log(acc)) : max;
	    
//		if (u>v)
//			return log1p(exp(v-u))+u;
//		else
//			return log1p(exp(u-v))+v;
	}

	
	public static double lbinom(int k, int n, double p) {
		return Binomial.density(k,n,p,true);
	}
	
	public static double ldtbbinom(int k, int n, double l, double u, double shape) {
		if (u<=l || k>n) return Double.NEGATIVE_INFINITY;
		
		double s1 = FastMath.exp(shape);
		double s2 = FastMath.exp(-shape);
		double cx = FastMath.min(k, n);
		return lchoose(n, k)+logdiff(IncompleteBeta.libeta(u,cx+s1,n-cx+s2),IncompleteBeta.libeta(l,cx+s1,n-cx+s2))-logdiff(IncompleteBeta.libeta(u,s1,s2),IncompleteBeta.libeta(l,s1,s2));
	}
	
	public static void rbinom(int n, double p, RandomEngine rnd, double[] re, int start, int end) {
		RandomState state = Binomial.create_random_state();
		for (int i = start; i<end; i++)
			re[i] = Binomial.random(n, p, rnd, state);
	}
	
	public static void rtbeta(double l, double u, double s1, double s2, RandomEngine rnd, double[] re, int start, int end) {
		// this is numerically super unstable (as qu tends to be close to 0)
//		double ql=Beta.cumulative(l,s1,s2,true,false);
//		double qu=Beta.cumulative(u,s1,s2,true,false);
//		for (int i=start; i<end; i++)
//			re[i] = Beta.quantile(Uniform.random(ql, qu, rnd), s1, s2, true, false);
		// do the direct inversion instead
		for (int i=start; i<end; i++)
			re[i] = Beta.quantile(
						lse( FastMath.log(Uniform.random(0, 1, rnd))+logdiff( IncompleteBeta.libeta(u, s1, s2), IncompleteBeta.libeta(l, s1, s2) ) , IncompleteBeta.libeta(l, s1, s2) )
						-MathFunctions.lbeta(s1,s2),
						s1,s2,true,true);
	}
	
	public static double rtbeta(double l, double u, double s1, double s2, RandomEngine rnd) {
		return Beta.quantile(
						lse( FastMath.log(Uniform.random(0, 1, rnd))+logdiff( IncompleteBeta.libeta(u, s1, s2), IncompleteBeta.libeta(l, s1, s2) ) , IncompleteBeta.libeta(l, s1, s2) )
						-MathFunctions.lbeta(s1,s2),
						s1,s2,true,true);
	}
	
	public static double rtbeta(double l, double u, double s1, double s2, double unif) {
		return Beta.quantile(
						lse( FastMath.log(unif)+logdiff( IncompleteBeta.libeta(u, s1, s2), IncompleteBeta.libeta(l, s1, s2) ) , IncompleteBeta.libeta(l, s1, s2) )
						-MathFunctions.lbeta(s1,s2),
						s1,s2,true,true);
	}

	public static void rtbeta2(double l, double u, double s1, double s2, RandomEngine rnd, double[] re, int start, int end) {
		for (int i=start; i<end; i++)
			re[i] = Beta.quantile(
						lse( FastMath.log(Uniform.random(0, 1, rnd))+logdiff( Beta.cumulative(u, s1, s2, true, true), Beta.cumulative(l, s1, s2, true, true) ) , Beta.cumulative(l, s1, s2, true, true) ),
						s1,s2,true,true);
	}

	
	public static void rtbbinom(int n, double l, double u, double shape, RandomEngine rnd, double[] re, int start, int end) {
		double s1=FastMath.exp(shape);
		double s2=FastMath.exp(-shape);
		rtbeta(l,u,s1,s2, rnd, re, start, end);
		RandomState state = Binomial.create_random_state();
		for (int i=start; i<end; i++)
			re[i] = Binomial.random(n, re[i], rnd, state);
	}

	
	public abstract double computeNewComponentExpectedProbability(MaximumLikelihoodParametrization par);
	public abstract double logLik(int k, int n, double[] param);
	public abstract double logLikOld(int k, int n, double[] param);
	public abstract double logLikNew(int k, int n, double[] param);
	
	@Override
	public double logLik(KNMatrixElement[] data, double[] param) {
		double re = 0;
		for (KNMatrixElement e : data)
			re+=logLik(e.k, e.n, param)*e.count;
		return re;
	}
	
	protected double[] sample(MaximumLikelihoodParametrization par, int count, int n, RandomEngine rnd) {
		int n1=(int) Binomial.random(count, par.getParameter(0), rnd);
		double[] re = new double[count];
		sampleFromComponent2(par,n,rnd, re,0,n1);
		sampleFromComponent1(par,n,rnd, re,n1,count);
		return re;
	}
	
	protected abstract void sampleFromComponent1(MaximumLikelihoodParametrization par, int n, RandomEngine rnd, double[] re, int start, int end);
	protected abstract void sampleFromComponent2(MaximumLikelihoodParametrization par, int n, RandomEngine rnd, double[] re, int start, int end);

	public KNMatrixElement[] sample(MaximumLikelihoodParametrization par, int[] nvec, RandomEngine rnd) {
		ArrayList<KNMatrixElement> re = new ArrayList<>();
		for (int n=1; n<nvec.length; n++) {
			KNMatrixElement[] tab = new KNMatrixElement[n+1];
			for (int i=0; i<tab.length; i++)
				tab[i] = new KNMatrixElement(i, n, 0);
			for (double k : sample(par,nvec[n],n,rnd))
				tab[(int)k].count++;
			for (int i=0; i<tab.length; i++)
				if (tab[i].count>0)
					re.add(tab[i]);
		}
		return re.toArray(new KNMatrixElement[0]);
	}
	
	
	public MaximumLikelihoodParametrization setNtrByFirstMoment(MaximumLikelihoodParametrization par, KNMatrixElement[] data) {
		// set expected NTR for given shape and pmconv
		double num = 0; double denom = 0;
		for (KNMatrixElement e : data){
			num+=e.count*e.k;
			denom+=e.count*e.n;
		}
		
		par.setParameter(0, Math.max(0.05, Math.min(0.9, num/denom/computeNewComponentExpectedProbability(par))));
		return par;
	}


	public static KNMatrixElement[] dataFromStatFile(String path, String condition, int subread, MetabolicLabelType label) throws IOException {
		return EI.lines(path).split('\t').skip(1)
			.filter(a->a[0].equals(condition) && a[1].equals(subread+"") && (a[2].equals(label.toString()) || a[2].equals(label.name())))
			.map(a->new KNMatrixElement(Integer.parseInt(a[3]), Integer.parseInt(a[4]), Double.parseDouble(a[5])))
			.toArray(KNMatrixElement.class);
	}
	
	public static HashMap<KnMatrixKey,KNMatrixElement[]> dataFromStatFile(String path) throws IOException {
		HashMap<KnMatrixKey,KNMatrixElement[]> re = new HashMap<>();
		
		for (String[][] block : EI.lines(path)
									.split('\t')
									.skip(1)
									.multiplexUnsorted((a,b)->a[0].equals(b[0]) && a[1].equals(b[1]) && a[2].equals(b[2]), String[].class)
									.loop()) {
			KnMatrixKey key = new KnMatrixKey(block[0][0],Integer.parseInt(block[0][1]), block[0][2]);
			KNMatrixElement[] value = EI.wrap(block)
					.map(a->new KNMatrixElement(Integer.parseInt(a[3]), Integer.parseInt(a[4]), Double.parseDouble(a[5])))
					.toArray(KNMatrixElement.class);
			re.put(key, value);
		}
			
		return re;
	}
	
//	public static void dataToStatFile(String path, KNMatrixElement[] elements, String condition, int subread, MetabolicLabelType label) throws IOException {
//		LineOrientedFile file = new LineOrientedFile(path);
//		boolean ex = file.exists();
//		LineWriter out = file.append();
//		if (!ex) 
//			out.writeLine("Condition\tSubread\tLabel\tk\tn\tCount");
//		EI.wrap(elements).map(e->String.format("%s\t%d\t%s\t%d\t%d\t%.0f", condition,subread,label.toString(),e.k,e.n,e.count)).print(out);
//		out.close();
//	}
	
	public static int[] nVectorFromStatFile(String path, String condition, int subread, MetabolicLabelType label) throws IOException {
		DoubleArrayList re = new DoubleArrayList();
		for (KNMatrixElement e : dataFromStatFile(path, condition, subread, label))
			re.increment(e.n, e.count);
		return re.toNumericArray().toIntArray();
	}


}
