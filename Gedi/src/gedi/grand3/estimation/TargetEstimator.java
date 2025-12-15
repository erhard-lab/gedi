package gedi.grand3.estimation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.BrentOptimizer;
import org.apache.commons.math3.optim.univariate.SearchInterval;
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction;
import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import gedi.grand3.estimation.models.Grand3BinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.SubreadKNKey;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.functions.EI;
import gedi.util.math.stat.inference.mixture.BiMixtureModelData;
import gedi.util.math.stat.inference.mixture.BiMixtureModelEstimator;
import gedi.util.math.stat.inference.mixture.BiMixtureModelResult;
import gedi.util.math.stat.inference.ml.NelderMeadMaximumLikelihoodEstimator;
import gedi.util.userInteraction.progress.Progress;

public class TargetEstimator {

	Grand3BinomialMixtureModel binomModel;
	Grand3TruncatedBetaBinomialMixtureModel tbbinomModel;
	
	
	NelderMeadMaximumLikelihoodEstimator estimator;
	
	// s,t,i (subread,type, sample)
	OldNewLogLikCache[][][] logLikComputerBinom;
	OldNewLogLikCache[][][] logLikComputerTbBinom;
	// shape,s,t,i (fixed shape parameter,subread,type, sample)
	OldNewLogLikCache[][][][] logLikComputerTbBinomShape;
	
	private ModelStructure[][][] models;
	private int[] conditionToSampleId;
	
	private int numSubreads = 0;
	private int numTypes = 0;
	private int numSamples = 0;
	private double[] shapes;
	private int[][] conditionMapping;
	
	private boolean fitMix;
	
	public TargetEstimator(ExperimentalDesign design, ModelStructure[][][] models, int[][] conditionMapping, int[] conditionToSampleId, double shapeStep, boolean fitMix) {
		this.models = models;
		this.conditionMapping = conditionMapping;
		this.conditionToSampleId = conditionToSampleId;
		this.fitMix = fitMix;
		for (ModelStructure[][] l2 : models)
			for (ModelStructure[] l1 : l2)
				for (ModelStructure m : l1) 
					if (m!=null) {
						numSubreads = Math.max(numSubreads, m.getSubread()+1);
						numTypes = Math.max(numTypes, m.getType()+1);
					}
		numSamples = design.getNumSamples();
		
		binomModel = new Grand3BinomialMixtureModel();
		tbbinomModel = new Grand3TruncatedBetaBinomialMixtureModel();
		
		estimator = new NelderMeadMaximumLikelihoodEstimator();
		estimator.setNthreads(0);
		
		double bottom = tbbinomModel.createPar().getLowerBounds("shape");
		double top = tbbinomModel.createPar().getUpperBounds("shape");
		shapes = new double[(int) (Math.ceil((top-bottom)/shapeStep))];
		for (int sh=0; sh<shapes.length; sh++) 
			shapes[sh] = Math.min(bottom+sh*shapeStep, top);
	}
	
	public void fillLikelihoodCache(Supplier<Progress> progress, int nthreads) {
		
		logLikComputerBinom = new OldNewLogLikCache[numSubreads][numTypes][numSamples];
		logLikComputerTbBinom = new OldNewLogLikCache[numSubreads][numTypes][numSamples];
		logLikComputerTbBinomShape = new OldNewLogLikCache[shapes.length][numSubreads][numTypes][numSamples];
		
		EI.seq(0, numSubreads)
			.unfold(s->EI.seq(0,numTypes).map(t->new int[]{s,t}))
			.unfold(p->EI.seq(0, numSamples).map(i->new int[]{p[0],p[1],i}))
			.unfold(p->EI.seq(0, 2+shapes.length).map(sh->new int[]{p[0],p[1],p[2],sh}))
			.progress(progress.get(),numSubreads*numTypes*numSamples*(2+shapes.length), ssti->"Caching likelihoods")
			.filter(sti->models[sti[0]][sti[1]][sti[2]]!=null)
			.parallelized(nthreads, 10, ei->ei.map(stish->{
				int s = stish[0];
				int t = stish[1];
				int i = stish[2];
				int sh = stish[3];
				if (sh==0)
					logLikComputerBinom[s][t][i] = new OldNewLogLikCache(binomModel,models[s][t][i].getBinom().getParameter());
				else if (sh==1)
					logLikComputerTbBinom[s][t][i] = new OldNewLogLikCache(tbbinomModel,models[s][t][i].getTBBinom().getParameter());
				else {
					sh = sh-2;
					double[] par = models[s][t][i].getTBBinom().getParameter().clone();
					par[models[s][t][i].getTBBinom().getParameterIndex("shape")] = shapes[sh];
					// ntr is not adapted, because not needed!
					logLikComputerTbBinomShape[sh][s][t][i] = new OldNewLogLikCache(tbbinomModel,par);
				}
				return 1;
			})).drain();
		
	}
	
	public int[][] getIndexMapping() {
		return conditionMapping;
	}
	
	public int getNumOutputConditions() {
		return conditionToSampleId.length;
	}
	
	public int getNumSamples() {
		return numSamples;
	}
	public int getNumSubreads() {
		return numSubreads;
	}
	public int getNumTypes() {
		return numTypes;
	}
	
	public boolean hasModel(int type, int conditionIndex) {
		for (int s = 0; s<logLikComputerBinom.length; s++)
			if (logLikComputerBinom[s][type][conditionToSampleId[conditionIndex]]!=null) 
				return true;
		return false;
	}

	public double getShape(int type, int conditionIndex) {
		for (int s = 0; s<logLikComputerBinom.length; s++)
			if (logLikComputerBinom[s][type][conditionToSampleId[conditionIndex]]!=null) 
				return models[s][type][conditionToSampleId[conditionIndex]].getTBBinom().getParameter("shape");
		throw new RuntimeException();
	}

	/**
	 * i is the condition index!
	 * @param s
	 * @param t
	 * @param i
	 * @return
	 */
	public OldNewLogLikCache getBinomialLogLikComputer(int s, int t, int i) {
		return logLikComputerBinom[s][t][conditionToSampleId[i]];
	}

	public OldNewLogLikCache getTbBinomialLogLikComputer(int s, int t, int i) {
		return logLikComputerTbBinom[s][t][conditionToSampleId[i]];
	}

	public OldNewLogLikCache getTbBinomialLogLikComputer(double shape, int s, int t, int i) {
		// find closest shapeindex
		int idx = Arrays.binarySearch(shapes, shape);
		if (idx==-1) idx = 0;
		else if (idx==-shapes.length-1) idx = shapes.length-1;
		else if (idx<0){
			double low = shapes[-idx-2];
			double high = shapes[-idx-1];
			idx = shape-low<high-shape?-idx-2:-idx-1;
		}
		return logLikComputerTbBinomShape[idx][s][t][conditionToSampleId[i]];
	}
	
	/**
	 * Returns the estimated shape, and sets the parameter into re
	 * @param t
	 * @param i
	 * @param llc
	 * @param re
	 * @return
	 */
	public double estimateTargetShape(int t, int i, HashMap<SubreadKNKey, AutoSparseDenseDoubleArrayCollector> llc,
			TargetEstimationResult re) {
		if (!hasModel(t,i)) throw new RuntimeException();
		
		ArrayList<SubreadKNKey> keys = new ArrayList<SubreadKNKey>();
		DoubleArrayList counts = new DoubleArrayList(); 
		for (SubreadKNKey key : llc.keySet()) {
			double cc = llc.get(key).get(i);
			if (cc>0) {
				keys.add(key);
				counts.add(cc);
			}
		}
		BiMixtureModelData[] buff = new BiMixtureModelData[keys.size()];
		for (int ii=0; ii<buff.length; ii++)
			buff[ii] = new BiMixtureModelData(0,0,0);
		
		double globLL = computeLL(getShape(t, i), keys, counts, buff, i, t);
		UnivariateFunction ufun = a->computeLL(a, keys, counts, buff, i, t);
		UnivariatePointValuePair opt = new BrentOptimizer(1E-2,1E-2).optimize(
				GoalType.MAXIMIZE,
				new SearchInterval(-2, 5),
				new UnivariateObjectiveFunction(ufun),
				new MaxEval(10000)
				);
		re.setTargetShapeInfo(t,i, opt.getPoint(),opt.getValue(),globLL);
		
		return opt.getPoint();
	}

	private double computeLL(double shape, ArrayList<SubreadKNKey> keys, DoubleArrayList counts, BiMixtureModelData[] data, int i, int t) {
		
		for (int ii=0; ii<data.length; ii++)
			updateBiMixModelData(data[ii],s->getTbBinomialLogLikComputer(shape,s,t,i), keys.get(ii), counts.getDouble(ii));
		
		BiMixtureModelEstimator model = new BiMixtureModelEstimator(data);
		return model.loglik(model.computeMAP());
	}
	
	private BiMixtureModelData updateBiMixModelData(BiMixtureModelData re, IntFunction<OldNewLogLikCache> llcomp, SubreadKNKey key, double count) {
		double oldLL = 0;
		double newLL = 0;
		for (int s=0; s<key.k.length; s++) {
			OldNewLogLikCache llc = llcomp.apply(s);
			if (llc!=null) {
				oldLL+=llc.getOldLogLik(key.k[s], key.n[s]);
				newLL+=llc.getNewLogLik(key.k[s], key.n[s]);
			}
		}
		return re.set(newLL, oldLL, count);
	}

	public void estimateMixtures(int t, int i, HashMap<SubreadKNKey, AutoSparseDenseDoubleArrayCollector> llc,
			TargetEstimationResult rec, double ci, boolean betaApprox, double shape) {
		if (!hasModel(t,i)) throw new RuntimeException();
		
		// collect all data
		ArrayList<BiMixtureModelData> binomMixData = new ArrayList<>();
		ArrayList<BiMixtureModelData> tbbinomMixData = new ArrayList<>();
		ArrayList<BiMixtureModelData> tbbinomShapeMixData = new ArrayList<>();
		for (SubreadKNKey key : llc.keySet()) {
			
			double count = llc.get(key).get(i);
			if (count>0) {
				// compute old and new ll for the three models
				binomMixData.add(createBiMixModelData(s->getBinomialLogLikComputer(s,t,i), key, count));
				tbbinomMixData.add(createBiMixModelData(s->getTbBinomialLogLikComputer(s,t,i), key, count));
				if (!Double.isNaN(shape))
					tbbinomShapeMixData.add(createBiMixModelData(s->getTbBinomialLogLikComputer(shape,s,t,i), key, count));
//				if (i==17)
//					System.out.println(key+" "+createBiMixModelData(s->getBinomialLogLikComputer(s,t,i), key, count));
			}
		}
		
		// estimate mixtures 
		if (!binomMixData.isEmpty()) {
			BiMixtureModelEstimator model = new BiMixtureModelEstimator(binomMixData.toArray(new BiMixtureModelData[binomMixData.size()]));
			model.setFitMix(fitMix);
			BiMixtureModelResult result = model.estimate(ci, betaApprox);
			rec.setTargetEstimateBinom(t,i,result);
			
			model = new BiMixtureModelEstimator(tbbinomMixData.toArray(new BiMixtureModelData[tbbinomMixData.size()]));
			model.setFitMix(fitMix);
			result = model.estimate(ci, betaApprox);
			rec.setTargetEstimateTbBinom(t,i, result);
			
			if (!Double.isNaN(shape)) {
				model = new BiMixtureModelEstimator(tbbinomShapeMixData.toArray(new BiMixtureModelData[tbbinomShapeMixData.size()]));
				model.setFitMix(fitMix);
				result = model.estimate(ci, betaApprox);
				rec.setTargetEstimateTbBinomShape(t,i, result);
			}
		}
	}
	
	private BiMixtureModelData createBiMixModelData(IntFunction<OldNewLogLikCache> llcomp, SubreadKNKey key, double count) {
		double oldLL = 0;
		double newLL = 0;
		for (int s=0; s<key.k.length; s++) {
			OldNewLogLikCache llc = llcomp.apply(s);
			if (llc!=null) { // could happen if there are no models for one of the subreads
				oldLL+=llc.getOldLogLik(key.k[s], key.n[s]);
				newLL+=llc.getNewLogLik(key.k[s], key.n[s]);
			}
		}
		return new BiMixtureModelData(newLL, oldLL, count);
	}


}

