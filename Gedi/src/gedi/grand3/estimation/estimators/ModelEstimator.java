package gedi.grand3.estimation.estimators;

import java.io.IOException;
import java.util.HashMap;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.DoubleFunction;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.logging.Logger;

import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.models.Grand3Model;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.KnMatrixKey;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModelPriorDecorator;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import gedi.util.math.stat.inference.ml.NelderMeadMaximumLikelihoodEstimator;
import jdistlib.Normal;

public abstract class ModelEstimator {

	protected ReentrantLock profLock = new ReentrantLock();
	protected Function<int[], KNMatrixElement[]> stiToData;
	protected double[][][][] pre_perr;
	protected String[] subreads;
	protected ExperimentalDesign design;
	protected boolean useGaussianWithoutno4sU;
	protected int maxIter;
	
	public ModelEstimator(ExperimentalDesign design, Function<int[], KNMatrixElement[]> stiToData, double[][][][] pre_perr, String[] subreads, boolean useGaussianWithoutno4sU) {
		this.design = design;
		this.stiToData = stiToData;
		this.pre_perr = pre_perr;
		this.subreads = subreads;
		this.useGaussianWithoutno4sU = useGaussianWithoutno4sU;
	}
	
	public ExtendedIterator<ModelStructure> estimate(int[] sti) {
		int s = sti[0];
		int t = sti[1];
		int i = sti[2];
		KNMatrixElement[] data = stiToData.apply(sti);
		if (data==null) return EI.empty();
		
		MaximumLikelihoodParametrization binom = fitBinom(data, pre_perr[s][t][i],t);
		MaximumLikelihoodParametrization tbbinom = fitTBBinom(data, pre_perr[s][t][i],t);

		if (getPrior(t).equals(LogDensityPrior.Epanechnikov)) {
			double m = (pre_perr[s][t][i][0]+pre_perr[s][t][i][1])/2;
			double r = (pre_perr[s][t][i][1]-pre_perr[s][t][i][0])/2;
			binom.setLowerBounds("p.err", m-r);
			binom.setUpperBounds("p.err", m+r);
			tbbinom.setLowerBounds("p.err", m-r);
			tbbinom.setUpperBounds("p.err", m+r);
		}
		
		return EI.wrap(new ModelStructure(s, t, i, binom, tbbinom));
	}
	
	public void setMaxIter(int maxIter) {
		this.maxIter = maxIter;
	}
	
	
	protected abstract MaximumLikelihoodParametrization fitBinom(KNMatrixElement[] data, double[] pre_perr, int t);
	protected abstract MaximumLikelihoodParametrization fitTBBinom(KNMatrixElement[] data, double[] pre_perr, int t);
	public LogDensityPrior getPrior(int t) {
		if (!useGaussianWithoutno4sU) return LogDensityPrior.Epanechnikov;
		boolean hasno4sU = design.getSamplesNotHaving(design.getTypes()[t]).length>0;
		return hasno4sU?LogDensityPrior.Epanechnikov:LogDensityPrior.Gaussian;
	}
	
	
	public String getPriorPerr() {
		StringBuilder sb = new StringBuilder();
		for (int t=0; t<design.getTypes().length; t++) {
			if (sb.length()>0) sb.append(", ");
			sb.append(getPrior(t).name()).append(" (").append(design.getTypes()[t].toString()).append(")");
		}
		return sb.toString();
	}
	
	public abstract ModelStructure computeProfile(Logger logger, LineWriter profOut, ModelStructure model, int nthreads, ExperimentalDesign design, String[] subreads);


	public static enum LogDensityPrior {
		Epanechnikov {
			public double logdens(double x, double mean, double rad) {
				double nd = (x-mean)/rad;
				nd = nd*nd;
				if (nd>=1) return Double.NEGATIVE_INFINITY;
//				return Math.log(0.75*(1-nd)/rad);
				return Math.log1p(-nd);
			}
		},
		Gaussian {
			public double logdens(double x, double mean, double rad) {
				return Normal.density(x, mean, rad/2, true);
			}
		};
		
		public abstract double logdens(double x, double mean, double rad);
	}
	
	
	public static void main(String[] args) throws IOException {
		HashMap<KnMatrixKey, KNMatrixElement[]> dataMap = Grand3Model.dataFromStatFile("/home/erhard/test/iter_117.conversion.knmatrix.tsv.gz");
		KNMatrixElement[] data = dataMap.values().iterator().next();
		
		Grand3TruncatedBetaBinomialMixtureModel tbbinomModel = new Grand3TruncatedBetaBinomialMixtureModel();
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		estimator.setMaxIter(10000);
		
		
		double m = 0.0002;
		double r = 0.00002;
		int perr = 1;
		ToDoubleFunction<double[]> perr_prior = a->LogDensityPrior.Epanechnikov.logdens(a[perr],m,r);
		
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> deco = new MaximumLikelihoodModelPriorDecorator<>(tbbinomModel,perr_prior);

		
		System.out.println(deco.logLik(data, new double[] {0.1118,0.0002374,0.1184,-1.57}));
		System.out.println(deco.logLik(data, new double[] {0.1718,0.0002374,0.1184,-1.57}));
		System.out.println(deco.logLik(data, new double[] {0.0591,0.0002,0.05,4.75}));
		
		MaximumLikelihoodParametrization par = tbbinomModel.createPar()
				.setParameter("shape", -1.6)
				.setParameter("p.err", 0.0002);
		System.out.println(estimator.fit(tbbinomModel, par,data));
		par = tbbinomModel.createPar()
				.setParameter("shape", -1.6)
				.setParameter("p.err", 0.0002);
		System.out.println(estimator.fit(deco, par,data));
		par = tbbinomModel.createPar()
				.setParameter("p.mconv", 0.1196)
				.setParameter("shape", -1.6)
				.setParameter("p.err", 0.0002);
		tbbinomModel.setNtrByFirstMoment(par,data);
		System.out.println("pre");
		System.out.println(par);
		System.out.println(estimator.fit(deco, par,data));
//		
//		
//		par = tbbinomModel.createPar()
//				.setParameter("shape", -1.6)
//				.setParameter("p.err", 0.0002);
//		DoubleFunction<MaximumLikelihoodParametrization> fit1 = s->{
//			MaximumLikelihoodParametrization par1= tbbinomModel.createPar()
//					.setParameter("shape", s)
//					.setParameter("p.err", m);
//			tbbinomModel.setNtrByFirstMoment(par1,data);
//			MaximumLikelihoodParametrization re = estimator.fit(deco,par1,data);
//			return re;
//		};
//		
//		System.out.println(fit1.apply(-1.6));
//		
		FullModelEstimator est = new FullModelEstimator(null, null, null, null);
		est.useGaussianWithoutno4sU = false;
		System.out.println(est.fitTBBinom(data,new double[] {0,0.0004},0));

	}
	
	
}
