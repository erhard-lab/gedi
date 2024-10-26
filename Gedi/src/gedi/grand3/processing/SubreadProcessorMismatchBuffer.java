package gedi.grand3.processing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.Grand3ReadClassified;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.ParallelizedState;

public class SubreadProcessorMismatchBuffer implements Grand3ReadClassified, ParallelizedState<SubreadProcessorMismatchBuffer> {
	
	private int[][] mm;
	
	public SubreadProcessorMismatchBuffer(int numSubreads) {
		mm = new int[numSubreads][16];
	}
	
	private CompatibilityCategory cat;
	private Collection<String> targets;
	
	private ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> read;
	int distinct;
	private ReadCountMode mode;
	
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(cat).append(" ").append(mode).append("\n");
		sb.append(ArrayUtils.matrixToString(mm));
		return sb.toString();
	}
	
	@Override
	public void classified(Collection<String> targets, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read,
			CompatibilityCategory cat, ReadCountMode mode, boolean sense) {
		this.targets = targets;
		this.read = (ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>) read;
		this.mode = mode;
		this.cat = cat;
		if (!sense) 
			throw new RuntimeException("Cannot be!");
	}
	
	public void increment(int s, char g) {
		int gi = SequenceUtils.inv_nucleotides[g];
		if (gi>=0 && gi<4)
			mm[s][gi*4+gi]++;
	}
	public void increment(int s, char g, char r) {
		int gi = SequenceUtils.inv_nucleotides[g];
		int ri = SequenceUtils.inv_nucleotides[r];
		if (gi>=0 && gi<4 && ri>=0 && ri<4)
			mm[s][gi*4+ri]++;
	}

	public void reset() {
		for (int[] mm : this.mm)
			Arrays.fill(mm,0);
	}
	
	public int getTotal(int s, MetabolicLabelType label) {
		if (getCategory().reversesLabel())
			return getTotal(s,label.getGenomicReverse());
		return getTotal(s,label.getGenomic());
	}
	
	public int getMismatches(int s, MetabolicLabelType label) {
		if (getCategory().reversesLabel())
			return getMismatches(s,label.getGenomicReverse(),label.getReadReverse());
		return getMismatches(s,label.getGenomic(),label.getRead());
	}
	
	public int getTotal(int s, char g) {
		int gi = SequenceUtils.inv_nucleotides[g];
		if (gi>=0 && gi<4)
			return mm[s][gi*4+gi];
		return 0;
	}
	
	public int getMismatches(int s, char g, char r) {
		int gi = SequenceUtils.inv_nucleotides[g];
		int ri = SequenceUtils.inv_nucleotides[r];
		if (gi>=0 && gi<4 && ri>=0 && ri<4)
			return mm[s][gi*4+ri];
		return 0;
	}


	
	public CompatibilityCategory getCategory() {
		return cat;
	}
	
	public AutoSparseDenseDoubleArrayCollector count(AutoSparseDenseDoubleArrayCollector re, int k, int[] reindex) {
		ReadCountMode mode = this.mode.transformCounts(k);
		if (!mode.equals(ReadCountMode.No))
			return read.getData().addCountsForDistinct(distinct, re, reindex, mode);
		return re;
	}
	
	public double[] count(double[] re, int k, int[] reindex) {
		ReadCountMode mode = this.mode.transformCounts(k);
		if (!mode.equals(ReadCountMode.No))
			return read.getData().addCountsForDistinct(distinct, re, reindex, mode);
		return re;
	}
	
	public AutoSparseDenseDoubleArrayCollector count(AutoSparseDenseDoubleArrayCollector re, int[] reindex) {
		return read.getData().addCountsForDistinct(distinct, re, reindex, mode);
	}
	public AutoSparseDenseDoubleArrayCollector count(AutoSparseDenseDoubleArrayCollector re, int[][] reindex) {
		return read.getData().addCountsForDistinct(distinct, re, reindex, mode);
	}
	
	public double[] count(double[] re, int[] reindex) {
		return read.getData().addCountsForDistinct(distinct, re, reindex, mode);
	}

	@Override
	public SubreadProcessorMismatchBuffer spawn(int index) {
		return new SubreadProcessorMismatchBuffer(mm.length);
	}

	@Override
	public void integrate(SubreadProcessorMismatchBuffer other) {
	}
	
	public ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> getRead() {
		return read;
	}
	
	public int getDistinct() {
		return distinct;
	}
	
	public Collection<String> getTargets() {
		return targets;
	}

	public ReadCountMode getMode() {
		return mode;
	}

	public boolean isSense() {
		return true;
	}

	public CompatibilityCategory getCompatibilityCategory() {
		return cat;
	}


	
}