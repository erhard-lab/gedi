package gedi.grand3.estimation.estimators;

import java.util.concurrent.locks.ReentrantLock;
import java.util.function.Function;
import java.util.logging.Logger;

import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import jdistlib.Normal;

public abstract class ModelEstimator {

	protected ReentrantLock profLock = new ReentrantLock();
	protected Function<int[], KNMatrixElement[]> stiToData;
	protected double[][][][] pre_perr;
	protected String[] subreads;
	protected ExperimentalDesign design;
	private boolean useGaussianWithoutno4sU;
	
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
	
	
}
