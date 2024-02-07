package gedi.grand3.estimation.estimators;

import java.io.IOException;
import java.util.function.Function;
import java.util.logging.Logger;

import cern.colt.bitvector.BitVector;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.models.Grand3BinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3MaskedBinomialModel;
import gedi.grand3.estimation.models.Grand3MaskedTruncatedBetaBinomialModel;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.ArrayUtils;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.StringLineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModelPriorDecorator;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import gedi.util.math.stat.inference.ml.NelderMeadMaximumLikelihoodEstimator;
import jdistlib.Binomial;

public class MaskErrorComponentModelEstimator extends ModelEstimator {

	

	private double maxErrComponentFraction = 0.01;
	
	public MaskErrorComponentModelEstimator(ExperimentalDesign design,Function<int[], KNMatrixElement[]> stiToData, double[][][][] pre_perr,
			String[] subreads) {
		super(design,stiToData, pre_perr, subreads,false);
	}

	public double getMaxErrComponentFraction() {
		return maxErrComponentFraction;
	}

	public void setMaxErrComponentFraction(double maxErrComponentFraction) {
		this.maxErrComponentFraction = maxErrComponentFraction;
	}

	public String getPriorPerr() {
		StringBuilder sb = new StringBuilder();
		for (int t=0; t<design.getTypes().length; t++) {
			sb.append("Epanechnikov (").append(design.getTypes()[t].name()).append(")");
		}
		return sb.toString();
	}

	private KNMatrixElement[] computeOnlyConversionComponent(KNMatrixElement[] data,double errorProb) {
		
		BitVector keep = new BitVector(data.length);
		
		int from=0;
		for (int to=1; to<=data.length; to++) {
			if (to==data.length || data[to].n!=data[from].n) {
				int n=data[from].n;
				
				double s=0;
				for (int i=from; i<to; i++) 
					s+=data[i].count;
				
				for (int i=from; i<to; i++) {
					int k = data[i].k;
					double expect = Binomial.density(k, n, errorProb, false)*s;
					keep.putQuick(i,expect<maxErrComponentFraction*data[i].count);
					if (i>from)
						keep.putQuick(i, keep.getQuick(i) || keep.getQuick(i-1));
				}
				
				from=to;
			}
		}
		return ArrayUtils.restrict(data, keep);
	}

	private static class Grand3NtrMixtureModel implements MaximumLikelihoodModel<KNMatrixElement[]>{

		private double pconv;
		private double shape;
		
		private Grand3NtrMixtureModel(double pconv) {
			this.pconv = pconv;
			this.shape = Double.NaN;
		}
		private Grand3NtrMixtureModel(double pmconv, double shape) {
			this.pconv = pmconv;
			this.shape = shape;
		}
		@Override
		public double logLik(KNMatrixElement[] data, double[] param) {
			double re = 0;
			if (Double.isNaN(shape))
				for (int i=0; i<data.length; i++) {
					re+=Grand3BinomialMixtureModel.ldbinommix(data[i].k, data[i].n, param[0], param[1], pconv)*data[i].count;
				}
			else
				for (int i=0; i<data.length; i++) {
					re+=Grand3TruncatedBetaBinomialMixtureModel.ldtbbinommix(data[i].k, data[i].n, param[0], param[1], pconv,shape)*data[i].count;
				}
			return re;
		}

		@Override
		public MaximumLikelihoodParametrization createPar() {
			return new MaximumLikelihoodParametrization(
					"binomNtr",
					new String[] {"ntr","p.err"},
					new double[] {0.2,4E-4},
					new String[] {"%.4f","%.4g"},
					new double[] {0,0},
					new double[] {1,0.1}
				);
		}
	}

	protected MaximumLikelihoodParametrization fitBinom(KNMatrixElement[] data, double[] pre_perr, int t) {
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;

		KNMatrixElement[] masked = computeOnlyConversionComponent(data, m);
		
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		Grand3MaskedBinomialModel model = new Grand3MaskedBinomialModel();
		MaximumLikelihoodParametrization re = estimator.fit(model,model.createPar(0.05), masked);

		Grand3NtrMixtureModel ntrModel = new Grand3NtrMixtureModel(re.getParameter(0));
		int perr = 1;
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> deco = new MaximumLikelihoodModelPriorDecorator<>(ntrModel,a->LogDensityPrior.Epanechnikov.logdens(a[perr],m,r));
		MaximumLikelihoodParametrization ntr = estimator.fit(deco,deco.createPar().setParameter(1, m), data);
		
		return new Grand3BinomialMixtureModel()
				.createPar(ntr.getParameter(0), ntr.getParameter(1), re.getParameter(0))
				.integrateOptimization(re)
				.integrateOptimization(ntr)
				.addConvergence(re.getConvergence())
				.addConvergence(ntr.getConvergence());
		
	}
	
	protected MaximumLikelihoodParametrization fitBinomJointSubread(KNMatrixElement[][] data, MaximumLikelihoodParametrization[] para, boolean ci, double[][] pre_perr, int t) {
		
		throw new RuntimeException();
	}

	protected MaximumLikelihoodParametrization fitTBBinom(KNMatrixElement[] data, double[] pre_perr, int t) {
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;

		KNMatrixElement[] masked = computeOnlyConversionComponent(data, m);
		
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		Grand3MaskedTruncatedBetaBinomialModel model = new Grand3MaskedTruncatedBetaBinomialModel(m);
		MaximumLikelihoodParametrization re = estimator.fit(model,model.createPar(0.05,4.5), masked);
		//estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(model, re, masked);
		
		Grand3NtrMixtureModel ntrModel = new Grand3NtrMixtureModel(re.getParameter(0),re.getParameter(1));
		int perr = 1;
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> deco = new MaximumLikelihoodModelPriorDecorator<>(ntrModel,a->LogDensityPrior.Epanechnikov.logdens(a[perr],m,r));

		MaximumLikelihoodParametrization ntr = estimator.fit(deco,deco.createPar().setParameter(1, m), data);

		return new Grand3TruncatedBetaBinomialMixtureModel().createPar(ntr.getParameter(0), ntr.getParameter(1), re.getParameter(0), re.getParameter(1)).integrateOptimization(re).integrateOptimization(ntr);
	}
	
	protected MaximumLikelihoodParametrization fitTBBinomJointSubread(KNMatrixElement[][] data, MaximumLikelihoodParametrization[] para, boolean ci, double[][] pre_perr) {
		throw new RuntimeException();
	}
	
	
	public ModelStructure computeProfile(Logger logger, LineWriter profOut, ModelStructure model, int nthreads, ExperimentalDesign design, String[] subreads) {

		String pref = String.format("MaskedError\t%s\t%s\t%s\t",design.getSampleNameForSampleIndex(model.getSample()),subreads[model.getSubread()],design.getTypes()[model.getType()].toString());
		
		double[] pre_perr = this.pre_perr[model.getSubread()][model.getType()][model.getSample()];
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;
		
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		
		// part 1: p.conv or p.mconv/shape
		Grand3MaskedBinomialModel binomModel = new Grand3MaskedBinomialModel(2);
		Grand3MaskedTruncatedBetaBinomialModel tbbinomModel = new Grand3MaskedTruncatedBetaBinomialModel(m,2,3);
		
		KNMatrixElement[] data = stiToData.apply(new int[]{model.getSubread(),model.getType(),model.getSample()});
		KNMatrixElement[] masked = computeOnlyConversionComponent(data, m);
		
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(binomModel,model.getBinom(),masked,"ntr","p.err");
		if (!model.getBinom().isConverged())
			logger.severe("Error computing binom CIs: "+pref+model.getBinom().getConvergence());
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(tbbinomModel,model.getTBBinom(),masked,"ntr","p.err");
		if (!model.getTBBinom().isConverged())
			logger.severe("Error computing tbbinom CIs: "+pref+model.getTBBinom().getConvergence());
		
		// part 2: ntr and p.err
		int perr = 1;
		Grand3NtrMixtureModel ntrModelBinom = new Grand3NtrMixtureModel(model.getBinom().getParameter("p.conv"));
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> decob = new MaximumLikelihoodModelPriorDecorator<>(ntrModelBinom,a->LogDensityPrior.Epanechnikov.logdens(a[perr],m,r));
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(decob,model.getBinom(),data,"p.conv");
		
		Grand3NtrMixtureModel ntrModelTbBinom = new Grand3NtrMixtureModel(model.getTBBinom().getParameter("p.mconv"),model.getTBBinom().getParameter("shape"));
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> decoTB = new MaximumLikelihoodModelPriorDecorator<>(ntrModelTbBinom,a->LogDensityPrior.Epanechnikov.logdens(a[perr],m,r));
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(decoTB,model.getTBBinom(),data,"p.mconv","shape");

		
		if (profOut!=null) {
			try {
				StringLineWriter sout = new StringLineWriter();
				estimator.writeProfileLikelihoods(tbbinomModel,model.getTBBinom(), masked, sout, pref,"ntr","p.err");
				if (!model.getTBBinom().isConverged())
					logger.severe("Error computing binom profiles: "+pref+model.getTBBinom().getConvergence());
				else {
					estimator.writeProfileLikelihoods(decoTB,model.getTBBinom(), data, sout, pref,"p.mconv","shape");
					if (!model.getTBBinom().isConverged())
						logger.severe("Error computing tbbinom profiles: "+pref+model.getTBBinom().getConvergence());
					else {
						profLock.lock();
						profOut.write(sout.toString());
						profLock.unlock();
					}
				}
			} catch (IOException e) {
			}
		}
		return model;
	}

}
