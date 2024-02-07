package gedi.lfc.localtest;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.mutable.MutableDouble;

public class LocalTestResult {


	private ImmutableReferenceGenomicRegion<String> gene;
	private double pval;
	
	private MemoryIntervalTreeStorage<MutableDouble> localPvalues;

	public LocalTestResult(ImmutableReferenceGenomicRegion<String> gene, double pval,
			MemoryIntervalTreeStorage<MutableDouble> localPvalues) {
		this.gene = gene;
		this.pval = pval;
		this.localPvalues = localPvalues;
	}

	
	public ImmutableReferenceGenomicRegion<String> getGene() {
		return gene;
	}
	
	
	public double getPvalue() {
		return pval;
	}
	
	public MemoryIntervalTreeStorage<MutableDouble> getLocalPvalues() {
		return localPvalues;
	}
	
}
