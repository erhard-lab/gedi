package gedi.core.data.reads.subreads;

import cern.colt.bitvector.BitVector;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.HasSubreads;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.SequenceUtils;

public class SingleEndToSubreadsConverter implements ReadByReadToSubreadsConverter<DefaultAlignedReadsData> {
	
	
	private boolean debug = false;
	private static final String[] semantic = {"Sense","Antisense"};
	
	@Override
	public String[] getSemantic() {
		return semantic;
	}

	
	
	@Override
	public ToSubreadsConverter setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	
	@Override
	public ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> convert(
			ImmutableReferenceGenomicRegion<? extends DefaultAlignedReadsData> read, boolean sense, MismatchReporter reporter) {

		
		HasSubreads subr = read.getData().asSubreads(sense);
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(read.getData().getNumConditions(), read.getData().hasNonzeroInformation());
		fac.start();
		
		AlignedReadsData re = read.getData();
		for (int d=0; d<re.getDistinctSequences(); d++) {
			fac.add(read.getData(), d,v->null,true);
			
			int[] sub = new int[subr.getNumSubreads(d)*2-1];
			for (int i=0; i<subr.getNumSubreads(d); i++) {
				sub[i*2] = subr.getSubreadId(d, i);
				if (i*2+1<sub.length)
					sub[i*2+1] = subr.getSubreadEnd(d, i, -1);
			}
			fac.setSubread(sub, subr.getGapPositions(d).toIntArray());
			BitVector retained = new BitVector(re.getVariationCount(d));
		
			for (int v=0; v<re.getVariationCount(d); v++) {
				if (re.isMismatch(d, v)) {
					fac.addMismatch(
							sense?re.getMismatchPos(d, v):(read.getRegion().getTotalLength()-1-re.getMismatchPos(d, v)), 
							sense^re.isVariationFromSecondRead(d, v)?re.getMismatchGenomic(d, v).charAt(0):SequenceUtils.getDnaComplement(re.getMismatchGenomic(d, v).charAt(0)), 
							sense^re.isVariationFromSecondRead(d, v)?re.getMismatchRead(d, v).charAt(0):SequenceUtils.getDnaComplement(re.getMismatchRead(d, v).charAt(0)), 
							false);
					if (reporter!=null) 
						retained.putQuick(v, true);
				}
				else if (re.isDeletion(d, v)) {
					fac.addDeletion(
							sense?re.getMismatchPos(d, v):(read.getRegion().getTotalLength()-re.getDeletionPos(d, v)-re.getDeletion(d, v).length()), 
							sense^re.isVariationFromSecondRead(d, v)?re.getDeletion(d, v):SequenceUtils.getDnaReverseComplement(re.getDeletion(d, v)), 
							false);
				}
				else if (re.isInsertion(d, v)) {
					fac.addInsertion(
							sense?re.getMismatchPos(d, v):(read.getRegion().getTotalLength()-1-re.getInsertionPos(d, v)), 
							sense^re.isVariationFromSecondRead(d, v)?re.getInsertion(d, v):SequenceUtils.getDnaReverseComplement(re.getInsertion(d, v)), 
							false);
				}
			}
			// report retained/corrected mismatches
			if (reporter!=null) {
				reporter.startDistinct(read, d, p->false);
				for (int v=0; v<re.getVariationCount(d); v++)
					if (re.isMismatch(d, v))
						reporter.reportMismatch(v, false, true);
			}
					
		}
		
		fac.makeDistinct();
		SubreadsAlignedReadsData data = fac.createSubread();
		
		ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> sread = new ImmutableReferenceGenomicRegion<>(sense?read.getReference():read.getReference().toOppositeStrand(), read.getRegion(), data);
		if (debug) {
			System.out.println("Convert:");
			System.out.println(read);
			System.out.println(sread);
			System.out.println();
		}
		return sread;
		
	}

	public static void main(String[] args) {
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
		fac.start();
		fac.newDistinctSequence();
		fac.addMismatch(10, 'C', 'T', false);
		ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read = ImmutableReferenceGenomicRegion.parse("1+:0-130",fac.create());
		System.out.println(read);
		System.out.println(new PairedEndToSubreadsConverter(false).convert(read, true, null));
		System.out.println(new PairedEndToSubreadsConverter(false).convert(read, false, null));
		
	}
	

}
