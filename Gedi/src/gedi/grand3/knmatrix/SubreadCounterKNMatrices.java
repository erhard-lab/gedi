package gedi.grand3.knmatrix;

import java.io.File;
import java.io.IOException;
import java.util.function.Predicate;

import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.processing.SubreadCounter;
import gedi.grand3.processing.SubreadProcessorMismatchBuffer;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public class SubreadCounterKNMatrices implements SubreadCounter<SubreadCounterKNMatrices> {

	
	private boolean debug = false;
	public SubreadCounterKNMatrices setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}

	@Override
	public SubreadCounterKNMatrices spawn(int index) {
		return new SubreadCounterKNMatrices(categories,counter[0][0].getNumConditions(),reindex,counter.length,labels).setDebug(debug);
	}

	@Override
	public void integrate(SubreadCounterKNMatrices other) {
		for (int c=0; c<counter.length; c++)
			for (int s=0; s<counter[c].length; s++)
				counter[c][s].integrate(other.counter[c][s]);
	}
	
	public KNMatrix getMatrix(int subread, int label) {
		return counter[subread][label];
	}
	
	
	// subread,label
	private KNMatrix[][] counter;
	private MetabolicLabelType[] labels;
	private Predicate<CompatibilityCategory> categories; 
	private int[] reindex;
	
	
	public SubreadCounterKNMatrices(Predicate<CompatibilityCategory> categories, int numCond, int[] reindex, int numSubreads, MetabolicLabelType[] labels) {
		this.categories = categories;
		this.labels = labels;
		counter = new KNMatrix[numSubreads][this.labels.length];
		for (int i=0; i<numSubreads; i++)
			for (int j=0; j<labels.length; j++)
				counter[i][j] = new KNMatrix(numCond);
		this.reindex = reindex;
	}

	@Override
	public void count(SubreadProcessorMismatchBuffer buffer) {
		int k,n;
		if (categories.test(buffer.getCategory())) {
		
			for (int s=0; s<counter.length; s++) {
				for (int l=0; l<labels.length; l++) {
					
					n = buffer.getTotal(s, labels[l]);
					k = buffer.getMismatches(s, labels[l]);
					
					if (k>n) 
						throw new RuntimeException("Cannot be: "+buffer.getRead());
					if (n>0) {
						buffer.count(counter[s][l].get(k, n),reindex);
						if (debug) {
							System.out.println("Binom: "+labels[l]+" s="+s+" k="+k+" n="+n);
						}
					}
				}
			}
		}
	}
	
	public void write(File out, ExperimentalDesign design) throws IOException {
		LineWriter wr = new LineOrientedFile(out.getPath()).write();
		wr.writeLine("Condition\tSubread\tLabel\tk\tn\tCount");
		for (int cond=0; cond<counter[0][0].getNumConditions(); cond++)
			for (int s=0; s<counter.length; s++)
				for (int l=0; l<labels.length; l++) {
					String prefix = String.format("%s\t%d\t%s\t", 
							design.getSampleNameForSampleIndex(cond),
							s,labels[l]);
					counter[s][l].write(prefix,cond,wr);
				}
		
		wr.close();
	}
	

}


