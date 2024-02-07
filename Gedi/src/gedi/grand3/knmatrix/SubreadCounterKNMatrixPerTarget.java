package gedi.grand3.knmatrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.function.Predicate;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.estimation.TargetEstimationResult;
import gedi.grand3.estimation.TargetEstimator;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.processing.SubreadCounter;
import gedi.grand3.processing.SubreadProcessorMismatchBuffer;
import gedi.grand3.processing.TargetCounter;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;

public class SubreadCounterKNMatrixPerTarget implements SubreadCounter<SubreadCounterKNMatrixPerTarget>, TargetCounter<SubreadCounterKNMatrixPerTarget,TargetEstimationResult> {

	
	private boolean debug = false;
	public SubreadCounterKNMatrixPerTarget setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	@Override
	public SubreadCounterKNMatrixPerTarget spawn(int index) {
		SubreadCounterKNMatrixPerTarget re = new SubreadCounterKNMatrixPerTarget(allCategories,useCat,targetEstimator,labels,subreadsToUse).setDebug(debug);
		re.ci = ci;
		re.betaApprox = betaApprox;
		return re;
	}

	@Override
	public void integrate(SubreadCounterKNMatrixPerTarget other) {
	}
	
	
	
	private double ci = 0.9;
	private boolean betaApprox = true;
	
	private ImmutableReferenceGenomicRegion<String> currentTarget; 
	private AutoSparseDenseDoubleArrayCollector[] total;
	// category,label
	private HashMap<SubreadKNKey,AutoSparseDenseDoubleArrayCollector>[][] llCounter;
	private SubreadKNKey llKey;
	
	private CompatibilityCategory[] allCategories;
	private Predicate<CompatibilityCategory> useCat;
	private TargetEstimator targetEstimator;
	private MetabolicLabelType[] labels;
	
	private boolean[] subreadsToUse;
	
	// TODO: Pool several conditions (e.g. cells,replicates) or genes? for estimating the shape parameter! 
	
	public SubreadCounterKNMatrixPerTarget(CompatibilityCategory[] allCategories, Predicate<CompatibilityCategory> useCat, TargetEstimator targetEstimator, MetabolicLabelType[] labels, boolean[] subreadsToUse) {
		this.allCategories = allCategories;
		this.useCat = useCat;
		this.targetEstimator = targetEstimator;
		this.labels = labels;
		this.subreadsToUse = subreadsToUse;
		int numCond = targetEstimator.getNumOutputConditions();
		total = new AutoSparseDenseDoubleArrayCollector[allCategories.length];

		llCounter = new HashMap[allCategories.length][targetEstimator.getNumTypes()];
		
		for (int c=0; c<allCategories.length; c++) 
			if (useCat.test(allCategories[c])) {
				total[c] = new AutoSparseDenseDoubleArrayCollector(numCond<50?numCond:10, numCond);
				for (int j=0; j<targetEstimator.getNumTypes(); j++) {
					llCounter[c][j] = new HashMap<>();
				}
			}
			else
				llCounter[c] = null;
		
		llKey = new SubreadKNKey(targetEstimator.getNumSubreads());
	}
	
	@Override
	public void startTarget(ImmutableReferenceGenomicRegion<String> target) {
		for (int a=0; a<llCounter.length; a++) 
			if (llCounter[a]!=null) {
				for (int b=0; b<llCounter[a].length; b++)
					llCounter[a][b].clear();
			}
		currentTarget = target;
		for (int a=0; a<total.length; a++) 
			if (total[a]!=null)
				total[a].clear();
	}

	@Override
	public void count(SubreadProcessorMismatchBuffer buffer) {
		if (!buffer.getTarget().equals(currentTarget))  throw new RuntimeException("Fatal error!");
		
		int numSub = targetEstimator.getNumSubreads();
		int numCond = targetEstimator.getNumOutputConditions();
		int k,n;
		if (useCat.test(buffer.getCategory())) {
			
			int c = buffer.getCategory().id();
			buffer.count(total[c], targetEstimator.getIndexMapping());
			
			for (int l=0; l<targetEstimator.getNumTypes(); l++) {
				
				for (int s=0; s<numSub; s++) {
							
					n = buffer.getTotal(s, labels[l]);
					k = buffer.getMismatches(s, labels[l]);
					if (subreadsToUse[s]) {
						llKey.k[s] = k;
						llKey.n[s] = n;
					} else {
						llKey.k[s] = 0;
						llKey.n[s] = 0;
					}
					
					if (k>n) 
						throw new RuntimeException("Cannot be: "+buffer.getRead());
					if (n>0) {
						if (debug) {
							System.out.println("Binom: "+labels[l]+" s="+s+" k="+k+" n="+n);
						}
					}
				}
				
				AutoSparseDenseDoubleArrayCollector co = llCounter[c][l].get(llKey);
				if (co==null) 
					llCounter[c][l].put(llKey.clone(), co=new AutoSparseDenseDoubleArrayCollector(numCond<50?numCond:10, numCond));
				buffer.count(co,targetEstimator.getIndexMapping());
				
			}
		}
	}
	

	@Override
	public List<TargetEstimationResult> getResultForCurrentTarget() {
		ArrayList<TargetEstimationResult> re = new ArrayList<TargetEstimationResult>();
		
		for (int c=0; c<llCounter.length; c++) {
			if (llCounter[c]==null) continue; // invalid category for target estimation
		
			TargetEstimationResult rec = new TargetEstimationResult(currentTarget,allCategories[c].getName(),labels.length);
			re.add(rec);
			
			int uc = c;
			total[c].iterate((i,count)->{
				rec.set(i, count);
				
				for (int t=0; t<labels.length; t++) {
					if (targetEstimator.hasModel(t, i)) {
						HashMap<SubreadKNKey, AutoSparseDenseDoubleArrayCollector> llc = llCounter[uc][t];
						// compute target specific shape 
						double shape = targetEstimator.estimateTargetShape(t, i, llc, rec);
						targetEstimator.estimateMixtures(t,i,llc,rec,ci,betaApprox,shape);
					}
				}
			});
			
			
		}
		return re;
		
	}
	
	
	

	

}


