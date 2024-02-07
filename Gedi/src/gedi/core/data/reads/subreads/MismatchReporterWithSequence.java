package gedi.core.data.reads.subreads;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.functions.IntToBooleanFunction;
import gedi.util.functions.ParallelizedState;

public abstract class MismatchReporterWithSequence<T extends MismatchReporterWithSequence<T>> implements MismatchReporter, ParallelizedState<T> {
	
	
	protected Genomic genomic;
	/**
	 * This is stateful and not thread safe. For parallel iteration, spawn copies per thread, and then collect them in the end!
	 * @param g
	 * @param targets 
	 */
	public MismatchReporterWithSequence(Genomic g) {
		this.genomic = g;
	}
	
	
	protected ImmutableReferenceGenomicRegion<String>  currentTarget;
	private ImmutableReferenceGenomicRegion<String>  currentTargetExtended;
	private char[] targetSeq;
	public void cacheRegion(ImmutableReferenceGenomicRegion<String> target, int extend) {
		this.currentTarget = target;
		this.currentTargetExtended = new ImmutableReferenceGenomicRegion<>(
				target.getReference(), 
				target.getRegion().extendAll(extend, extend).intersect(0, genomic.getLength(target.getReference().getName())),
				target.getData());
		
		targetSeq = genomic.getSequence(currentTargetExtended).toString().toUpperCase().toCharArray();
	}
	
	protected ImmutableReferenceGenomicRegion<? extends AlignedReadsData> currentRead;
	protected int currentDistinct;
	protected char[] readSeq;
	
	
	@Override
	public void startDistinct(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read,
			int distinct,IntToBooleanFunction overlapperGetter) {
		
		if (currentRead==null || currentRead.compareTo((ImmutableReferenceGenomicRegion)read)!=0) {
			// obtain sequence
			if (currentTargetExtended.getRegion().contains(read.getRegion())) {
				readSeq = SequenceUtils.extractSequence(currentTargetExtended.induce(read.getRegion()), targetSeq);
				if (!read.getReference().getStrand().equals(currentTargetExtended.getReference().getStrand()))
					SequenceUtils.getDnaReverseComplementInplace(readSeq);
			}
			else {
				readSeq = genomic.getSequence(read).toString().toUpperCase().toCharArray();
			}
		}
		
		currentRead = read;
		currentDistinct = distinct;
	}
	
	
	/**
	 * Should usually return true. If the variation in the read is incorrect, it will either return false or throw an Exception. 
	 * It will only be false, if the position in the genome is anything else than ACGT, and in the read it is N! (is the case for ERCC_138!)
	 *  
	 * @param variation
	 * @param inOverlap
	 * @param retained
	 * @return
	 */
	protected boolean checkSequence(int variation) {
		// test genomic sequence
		char mm = currentRead.getData().getMismatchGenomic(currentDistinct, variation).charAt(0);
		char genome = (currentRead.getData().isVariationFromSecondRead(currentDistinct, variation)?SequenceUtils.getDnaComplement(readSeq[currentRead.getData().getMismatchPos(currentDistinct, variation)]):readSeq[currentRead.getData().getMismatchPos(currentDistinct, variation)]);
		if (mm=='N' && genome!='A'&&genome!='C'&&genome!='G'&&genome!='T')
			return false;
		if (mm!=genome)
			throw new RuntimeException("Sequences do not match! This is a sign that references for read mapping and Grand3 are different!\n"+currentRead+"\nReference sequence: "+String.valueOf(readSeq)+"\n"+currentRead.getData().getVariation(currentDistinct, variation));
		return true;
	}
	

	

}
