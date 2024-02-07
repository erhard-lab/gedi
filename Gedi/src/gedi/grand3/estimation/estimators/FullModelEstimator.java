package gedi.grand3.estimation.estimators;

import java.io.IOException;
import java.util.Stack;
import java.util.function.DoubleFunction;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.logging.Logger;

import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.models.Grand3BinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3SubreadsJointBinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3SubreadsJointTruncatedBetaBinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.StringLineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModelPriorDecorator;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import gedi.util.math.stat.inference.ml.NelderMeadMaximumLikelihoodEstimator;

public class FullModelEstimator extends ModelEstimator {

	

	private double minShapeStep = 4;
	private double maxShapeStep = 0.5;
	
	
	public FullModelEstimator(ExperimentalDesign design, Function<int[], KNMatrixElement[]> stiToData, double[][][][] pre_perr,
			String[] subreads) {
		super(design,stiToData, pre_perr, subreads,true);
	}

	public void setMaxShapeStep(double maxShapeStep) {
		this.maxShapeStep = maxShapeStep;
	}
	
	public void setMinShapeStep(double minShapeStep) {
		this.minShapeStep = minShapeStep;
	}
	
	public MaximumLikelihoodParametrization[] estimateJoint(ModelStructure[] block, boolean profile) {
		KNMatrixElement[][] data = new KNMatrixElement[block.length][];
		MaximumLikelihoodParametrization[] binomPara = new MaximumLikelihoodParametrization[block.length];
		MaximumLikelihoodParametrization[] tbbinomPara = new MaximumLikelihoodParametrization[block.length];
		double[][] tpre_perr = new double[block.length][];
		
		int[] ts = new int[block.length];
		
		for (int j=0; j<block.length; j++) {
			data[j] = stiToData.apply(new int[] {block[j].getSubread(),block[j].getType(),block[j].getSample()});
			binomPara[j] = block[j].getBinom();
			tbbinomPara[j] = block[j].getTBBinom();
			tpre_perr[j] = pre_perr[block[j].getSubread()][block[j].getType()][block[j].getSample()];
			ts[j] = block[j].getType();
		}
		
		
		MaximumLikelihoodParametrization model = fitBinomJointSubread(data, binomPara, profile,tpre_perr,ts);
		MaximumLikelihoodParametrization tbmodel = fitTBBinomJointSubread(data, tbbinomPara, profile,tpre_perr,ts);
		return new MaximumLikelihoodParametrization[] {model,tbmodel};
	}
	
	protected MaximumLikelihoodParametrization fitBinom(KNMatrixElement[] data, double[] pre_perr, int t) {
		Grand3BinomialMixtureModel binomModel = new Grand3BinomialMixtureModel();
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;
		int perr = 1;
		ToDoubleFunction<double[]> perr_prior = a->getPrior(t).logdens(a[perr],m,r);
		

		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> deco = new MaximumLikelihoodModelPriorDecorator<>(binomModel,perr_prior);

		MaximumLikelihoodParametrization par = binomModel.createPar()
				.setParameter("p.err", m);

		binomModel.setNtrByFirstMoment(par,data);
		MaximumLikelihoodParametrization re = estimator.fit(deco,par, data);
//		re.setLowerBounds("p.err", m-r);
//		re.setUpperBounds("p.err", m+r); // not sensible for gaussian prior!
		return re;
		
//		MaximumLikelihoodParametrization par = binomModel.createPar()
//				.setParameter("p.err", (pre_perr[0]+pre_perr[1])/2);
//		binomModel.setNtrByFirstMoment(par,data);
//		return estimator.fit(binomModel,par, data, "p.err");
	}
	
	protected MaximumLikelihoodParametrization fitBinomJointSubread(KNMatrixElement[][] data, MaximumLikelihoodParametrization[] para, boolean ci, double[][] pre_perr, int[] t) {
		
		Grand3BinomialMixtureModel binomModel = new Grand3BinomialMixtureModel();
		
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]>[] deco = new MaximumLikelihoodModelPriorDecorator[para.length];
		for (int i=0; i<para.length; i++) {
			double m = (pre_perr[i][0]+pre_perr[i][1])/2;
			double r = (pre_perr[i][1]-pre_perr[i][0])/2;
			
			int perr = para[i].getParameterIndex("p.err");
			int ii=i;
			ToDoubleFunction<double[]> perr_prior = a->getPrior(t[ii]).logdens(a[perr],m,r);
			
			deco[i] = new MaximumLikelihoodModelPriorDecorator<>(binomModel,perr_prior);
		}
		
		Grand3SubreadsJointBinomialMixtureModel model = new Grand3SubreadsJointBinomialMixtureModel(deco);
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();

		
		MaximumLikelihoodParametrization par = model.createPar(para);
		MaximumLikelihoodParametrization re = estimator.fit(model,par, data);
		if (ci) estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(model, re, data);
		
		model.distributeTo(par,para);
		return re;
	}

	protected MaximumLikelihoodParametrization fitTBBinom(KNMatrixElement[] data, double[] pre_perr, int t) {
		
		Grand3TruncatedBetaBinomialMixtureModel tbbinomModel = new Grand3TruncatedBetaBinomialMixtureModel();
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;
		int perr = 1;
		ToDoubleFunction<double[]> perr_prior = a->getPrior(t).logdens(a[perr],m,r);
		
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> deco = new MaximumLikelihoodModelPriorDecorator<>(tbbinomModel,perr_prior);

		DoubleFunction<MaximumLikelihoodParametrization> fit1 = s->{
			MaximumLikelihoodParametrization par = tbbinomModel.createPar()
					.setParameter("shape", s)
					.setParameter("p.err", m);
			tbbinomModel.setNtrByFirstMoment(par,data);
			MaximumLikelihoodParametrization re = estimator.fit(deco,par,data);
			return re;
		};
		
		double bottom = tbbinomModel.createPar().getLowerBounds("shape")+maxShapeStep/2;
		double top = tbbinomModel.createPar().getUpperBounds("shape")-maxShapeStep/2;
		
		Stack<ShapeBracket> toFit = new Stack<ShapeBracket>(); 
		toFit.add(new ShapeBracket(bottom,top,fit1.apply(bottom),fit1.apply(top)));
		
		MaximumLikelihoodParametrization bestModel = toFit.peek().getBetter();
		while (!toFit.isEmpty()) {
			ShapeBracket bracket = toFit.pop();
			MaximumLikelihoodParametrization better = bracket.getBetter();
			MaximumLikelihoodParametrization other = bracket.getWorse();
			
			if ((bracket.topShape-bracket.bottomShape>minShapeStep) || (bracket.topShape-bracket.bottomShape>maxShapeStep && !better.isSimilar(other, 1E-3))) {
				double centerShape = (bracket.topShape+bracket.bottomShape)/2;
				MaximumLikelihoodParametrization center = fit1.apply(centerShape);
				if (center.getLogLik()>bestModel.getLogLik())
					bestModel = center;
				toFit.add(new ShapeBracket(bracket.bottomShape, centerShape, bracket.bottomModel, center));
				toFit.add(new ShapeBracket(centerShape, bracket.topShape, center, bracket.topModel));
			}
		}
		MaximumLikelihoodParametrization tbbinom = bestModel;
//		tbbinom.setLowerBounds("p.err", m-r);
//		tbbinom.setUpperBounds("p.err", m+r); // not sensible for gaussian prior!
		
		return tbbinom;
	}
	
	protected MaximumLikelihoodParametrization fitTBBinomJointSubread(KNMatrixElement[][] data, MaximumLikelihoodParametrization[] para, boolean ci, double[][] pre_perr, int[] t) {
		
		Grand3TruncatedBetaBinomialMixtureModel tbbinomModel = new Grand3TruncatedBetaBinomialMixtureModel();
		
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]>[] deco = new MaximumLikelihoodModelPriorDecorator[para.length];
		for (int i=0; i<para.length; i++) {
			double m = (pre_perr[i][0]+pre_perr[i][1])/2;
			double r = (pre_perr[i][1]-pre_perr[i][0])/2;
			
			int perr = para[i].getParameterIndex("p.err");
			int ii=i;
			ToDoubleFunction<double[]> perr_prior = a->getPrior(t[ii]).logdens(a[perr],m,r);
						
			deco[i] = new MaximumLikelihoodModelPriorDecorator<>(tbbinomModel,perr_prior);
		}
		
		Grand3SubreadsJointTruncatedBetaBinomialMixtureModel model = new Grand3SubreadsJointTruncatedBetaBinomialMixtureModel(deco);
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		
		MaximumLikelihoodParametrization par = model.createPar(para);
		MaximumLikelihoodParametrization re = estimator.fit(model,par, data);
		if (ci) estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(model, re, data);
		
		model.distributeTo(par,para);
		return re;
	}
	
	
	public ModelStructure computeProfile(Logger logger, LineWriter profOut, ModelStructure model, int nthreads, ExperimentalDesign design, String[] subreads) {

		Grand3BinomialMixtureModel binomModel = new Grand3BinomialMixtureModel();
		Grand3TruncatedBetaBinomialMixtureModel tbbinomModel = new Grand3TruncatedBetaBinomialMixtureModel();
		NelderMeadMaximumLikelihoodEstimator estimator = new NelderMeadMaximumLikelihoodEstimator();
		
		
		double[] pre_perr = this.pre_perr[model.getSubread()][model.getType()][model.getSample()];
		
		
		double m = (pre_perr[0]+pre_perr[1])/2;
		double r = (pre_perr[1]-pre_perr[0])/2;
		int perr = 1;

		int t = model.getType();
		ToDoubleFunction<double[]> perr_prior = a->getPrior(t).logdens(a[perr],m,r);
			
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> bdeco = new MaximumLikelihoodModelPriorDecorator<>(binomModel,perr_prior);
		MaximumLikelihoodModelPriorDecorator<KNMatrixElement[]> tbbdeco = new MaximumLikelihoodModelPriorDecorator<>(tbbinomModel,perr_prior);

		KNMatrixElement[] data = stiToData.apply(new int[]{model.getSubread(),model.getType(),model.getSample()});
		
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(bdeco,model.getBinom(),data);
		estimator.computeProfileLikelihoodPointwiseConfidenceIntervals(tbbdeco,model.getTBBinom(),data);
		if (profOut!=null) {
			try {
				String pref = String.format("Separate\t%s\t%s\t%s\t",design.getSampleNameForSampleIndex(model.getSample()),subreads[model.getSubread()],design.getTypes()[model.getType()].toString());
				StringLineWriter sout = new StringLineWriter();
				estimator.writeProfileLikelihoods(tbbdeco,model.getTBBinom(), data, sout, pref);
				if (!model.getTBBinom().isConverged())
					logger.severe(pref+model.getTBBinom().getConvergence());
				profLock.lock();
				profOut.write(sout.toString());
				profLock.unlock();
			} catch (IOException e) {
			}
		}
		return model;
	}

	private static class ShapeBracket {
		double bottomShape;
		double topShape;
		MaximumLikelihoodParametrization bottomModel;
		MaximumLikelihoodParametrization topModel;
		public ShapeBracket(double bottomShape, double topShape, MaximumLikelihoodParametrization bottomModel,
				MaximumLikelihoodParametrization topModel) {
			this.bottomShape = bottomShape;
			this.topShape = topShape;
			this.bottomModel = bottomModel;
			this.topModel = topModel;
		}
		public MaximumLikelihoodParametrization getBetter() {
			return bottomModel.getLogLik()>topModel.getLogLik()?bottomModel:topModel;
		}
		public MaximumLikelihoodParametrization getWorse() {
			return bottomModel.getLogLik()>topModel.getLogLik()?topModel:bottomModel;
		}
		
		
		
	}
}
