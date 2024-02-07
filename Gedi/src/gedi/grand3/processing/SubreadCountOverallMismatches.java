package gedi.grand3.processing;

import java.io.File;
import java.io.IOException;
import java.util.function.IntFunction;

import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.TargetCollection;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public class SubreadCountOverallMismatches implements SubreadCounter<SubreadCountOverallMismatches> {

	@Override
	public SubreadCountOverallMismatches spawn(int index) {
		return new SubreadCountOverallMismatches(counter.length, counter[0].length,reindex);
	}

	@Override
	public void integrate(SubreadCountOverallMismatches other) {
		for (int c=0; c<counter.length; c++)
			for (int s=0; s<counter[c].length; s++)
				for (int i=0; i<counter[c][s].length; i++)
					ArrayUtils.add(counter[c][s][i], other.counter[c][s][i]);
	}
	
	
	// cat,subread,mismatch,condition
	private double[][][][] counter;
	private int[] reindex;
	
	public SubreadCountOverallMismatches(int numCategories, int numSubreads, int[] reindex) {
		counter = new double[numCategories][numSubreads][16][ArrayUtils.max(reindex)+1];
		this.reindex = reindex;
	}

	@Override
	public void count(SubreadProcessorMismatchBuffer buffer) {
		for (int s=0; s<counter[0].length; s++) {
			for (char g : SequenceUtils.valid_nucleotides) {
				for (char r : SequenceUtils.valid_nucleotides){
					int k = buffer.getMismatches(s, g, r); // if g==r this is no of matches!
					if (k>0) {
						int idx = SequenceUtils.inv_nucleotides[g]*4+SequenceUtils.inv_nucleotides[r];
						buffer.count(counter[buffer.getCategory().id()][s][idx],k,reindex);
					}
				}
			}
		}
	}
	
	public void write(File out, IntFunction<String> nameFunction, IntFunction<String> condNameFunction, TargetCollection targets) throws IOException {
		CompatibilityCategory[] cats = targets.getCategories();
		
		LineWriter wr = new LineOrientedFile(out.getPath()).write();
		wr.writeLine("Category\tName\tCondition\tSubread\tGenomic\tRead\tFrequency");
		for (int cond=0; cond<counter[0][0][0].length; cond++)
			for (int c=0; c<counter.length; c++)
				for (int s=0; s<counter[c].length; s++)
					for (int g=0; g<4; g++)
						for (int r=0; r<4; r++)
							if (r!=g && counter[c][s][g*4+g][cond]>0)
								wr.writef("%s\t%s\t%s\t%d\t%s\t%s\t%.3g\n", 
										cats[c],nameFunction.apply(cond),
										condNameFunction.apply(cond),
										s,SequenceUtils.valid_nucleotides[g],SequenceUtils.valid_nucleotides[r],
										counter[c][s][g*4+r][cond]/counter[c][s][g*4+g][cond]
												);
		
		wr.close();
	}
	

	public void writePerr(File out, IntFunction<String> condNameFunction, TargetCollection targets) throws IOException {
		CompatibilityCategory[] cats = targets.getCategories();
		
		double[][][] counter = new double[this.counter[0][0][0].length][this.counter[0].length][16];

		
		for (int cond=0; cond<this.counter[0][0][0].length; cond++)
			for (int c=0; c<this.counter.length; c++)
				if (cats[c].useToEstimateGlobalParameters())
					for (int s=0; s<this.counter[c].length; s++)
						for (int gr=0; gr<16; gr++)
								counter[cond][s][gr]+=this.counter[c][s][gr][cond];
		
		LineWriter wr = new LineOrientedFile(out.getPath()).write();
		wr.writeLine("Condition\tSubread\tGenomic\tRead\tFrequency");
		for (int cond=0; cond<counter.length; cond++)
				for (int s=0; s<counter[cond].length; s++)
					for (int g=0; g<4; g++)
						for (int r=0; r<4; r++)
							if (r!=g && counter[cond][s][g*4+g]>0)
								wr.writef("%s\t%d\t%s\t%s\t%.3g\n", 
										condNameFunction.apply(cond),
										s,SequenceUtils.valid_nucleotides[g],SequenceUtils.valid_nucleotides[r],
										counter[cond][s][g*4+r]/counter[cond][s][g*4+g]
												);
		
		wr.close();
	}
	
}

