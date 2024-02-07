package gedi.grand3.processing;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.Grand3ReadClassified;
import gedi.util.functions.ParallelizedState;
import jdistlib.rng.MersenneTwister;
import jdistlib.rng.RandomEngine;

public class ResimulateState implements Grand3ReadClassified, ParallelizedState<ResimulateState> {

	
	private long seed;
	
	public ResimulateState(long seed) {
		this.seed = seed;
	}

	private String target;
	private CompatibilityCategory cat;
	
	
	@Override
	public void classified(ImmutableReferenceGenomicRegion<String> target, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read,
			CompatibilityCategory cat, ReadCountMode mode, boolean sense) {
		if (mode.equals(ReadCountMode.No) || cat==null || !cat.useToEstimateTargetParameters())
			this.target = null;
		else
			this.target = target.getData();
		this.cat = cat;
		
		if (!sense) 
			throw new RuntimeException("Cannot be!");
	}
	
	private ImmutableReferenceGenomicRegion<String> region;
	private char[] sequence;
	private RandomEngine rnd;
	public void setSequence(ImmutableReferenceGenomicRegion<String> region, char[] sequence) {
		this.region = region;
		this.sequence = sequence;
		rnd = new MersenneTwister(seed+region.getRegion().getStart());
	}

	public RandomEngine getRnd() {
		return rnd;
	}

	@Override
	public ResimulateState spawn(int index) {
		return new ResimulateState(seed);
	}

	@Override
	public void integrate(ResimulateState other) {
	}
	
	public String getTarget() {
		return target;
	}

	public CompatibilityCategory getCategory() {
		return cat;
	}

	public char[] getSequence() {
		return sequence;
	}

	public ImmutableReferenceGenomicRegion<String> getSequenceRegion() {
		return region;
	}

	
}