package gedi.core.data.reads.subreads;

import cern.colt.bitvector.BitVector;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.HasSubreads;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.functions.BiIntConsumer;

public class PairedEndToSubreadsConverter implements ReadByReadToSubreadsConverter<DefaultAlignedReadsData> {
	
	private boolean debug = false;
	private boolean alldeletions;
	
	private static final String[] semantic = {"Sense","Overlap","Antisense"};
	
	@Override
	public String[] getSemantic() {
		return semantic;
	}
	
	/**
	 * alldeletions means to keep also deletion only observed in one read. This makes sense for counting consistent mismatches 
	 * (as for grand3), as these bases should not be considered (there cant be a consistent mismatch in them, and they have not been
	 * sequenced, but are part of the mapped region).
	 * 
	 * @param alldeletions
	 * @param debug
	 * @return
	 */
	public PairedEndToSubreadsConverter(boolean alldeletions) {
		this.alldeletions = alldeletions;
	}
	
	public PairedEndToSubreadsConverter setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	
	private static int compareVariation(AlignedReadsData ard, int d, int v1, int v2) {
		int re = Integer.compare(ard.getMismatchPos(d, v1), ard.getMismatchPos(d, v2));
		if (re!=0) return re;
		if (ard.isMismatch(d,v1)!=ard.isMismatch(d,v2)) return ard.isMismatch(d,v1)?-1:1;
		if (ard.isInsertion(d,v1)!=ard.isInsertion(d,v2)) return ard.isInsertion(d,v1)?-1:1;
		if (ard.isDeletion(d,v1)!=ard.isDeletion(d,v2)) return ard.isDeletion(d,v1)?-1:1;
		if (ard.isSoftclip(d,v1)!=ard.isSoftclip(d,v2)) return ard.isSoftclip(d,v1)?-1:1;
		if (ard.isVariationFromSecondRead(d,v1)!=ard.isVariationFromSecondRead(d,v2)) return ard.isVariationFromSecondRead(d,v1)?1:-1;
		return 0;
	}
	
	@Override
	public ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> convert(
			ImmutableReferenceGenomicRegion<? extends DefaultAlignedReadsData> read, boolean sense, MismatchReporter reporter,
			BiIntConsumer usedTotal) {

		
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
				if (v>0 && compareVariation(re, d, v-1, v)>=0)
					throw new RuntimeException("Variations are not sorted: "+read);
				if (re.isMismatch(d, v)) {
					if ((! re.isPositionInOverlap(d, re.getMismatchPos(d, v)) ||
							(v>0 &&
							re.isMismatch(d, v-1) && 
							re.getMismatchPos(d, v-1)==re.getMismatchPos(d, v) &&
							re.getMismatchGenomic(d, v-1).charAt(0)==SequenceUtils.getDnaComplement(re.getMismatchGenomic(d, v).charAt(0)) &&
							re.getMismatchRead(d, v-1).charAt(0)==SequenceUtils.getDnaComplement(re.getMismatchRead(d, v).charAt(0)) &&
							re.isVariationFromSecondRead(d, v-1)!=re.isVariationFromSecondRead(d, v)
							))) {
						fac.addMismatch(
								sense?re.getMismatchPos(d, v):(read.getRegion().getTotalLength()-1-re.getMismatchPos(d, v)), 
								sense^re.isVariationFromSecondRead(d, v)?re.getMismatchGenomic(d, v).charAt(0):SequenceUtils.getDnaComplement(re.getMismatchGenomic(d, v).charAt(0)), 
								sense^re.isVariationFromSecondRead(d, v)?re.getMismatchRead(d, v).charAt(0):SequenceUtils.getDnaComplement(re.getMismatchRead(d, v).charAt(0)), 
								false);
						if (reporter!=null) {
							if (re.isPositionInOverlap(d, re.getMismatchPos(d, v)))
								retained.putQuick(v-1, true);
							retained.putQuick(v, true);
						}
					} 
				}
				else if (re.isDeletion(d, v)) {
					boolean secondOfPairedDeletion = (v>0 &&
							re.isDeletion(d, v-1) && 
							re.getDeletionPos(d, v-1)==re.getDeletionPos(d, v) &&
							re.getDeletion(d, v-1).equals(SequenceUtils.getDnaReverseComplement(re.getDeletion(d, v))) &&
							re.isVariationFromSecondRead(d, v-1)!=re.isVariationFromSecondRead(d, v)
							);
					if (!re.isPositionInOverlap(d, re.getDeletionPos(d, v)) || secondOfPairedDeletion || (alldeletions && !secondOfPairedDeletion)) {
						fac.addDeletion(
								sense?re.getDeletionPos(d, v):(read.getRegion().getTotalLength()-re.getDeletionPos(d, v)-re.getDeletion(d, v).length()), 
								sense^re.isVariationFromSecondRead(d, v)?re.getDeletion(d, v):SequenceUtils.getDnaReverseComplement(re.getDeletion(d, v)), 
								false);
					}
				}
				else if (re.isInsertion(d, v)) {
					if ((! re.isPositionInOverlap(d, re.getInsertionPos(d, v)) ||
							(v>0 &&
							re.isInsertion(d, v-1)  && 
							re.getInsertionPos(d, v-1)==re.getInsertionPos(d, v) &&
							re.getInsertion(d, v-1).equals(SequenceUtils.getDnaReverseComplement(re.getInsertion(d, v))) &&
							re.isVariationFromSecondRead(d, v-1)!=re.isVariationFromSecondRead(d, v)
							))) {
						fac.addInsertion(
								sense?re.getInsertionPos(d, v):(read.getRegion().getTotalLength()-1-re.getInsertionPos(d, v)), 
								sense^re.isVariationFromSecondRead(d, v)?re.getInsertion(d, v):SequenceUtils.getDnaReverseComplement(re.getInsertion(d, v)), 
								false);
					}
				}
			}
			// report retained/corrected mismatches
			if (reporter!=null) {
				int ud = d;
				reporter.startDistinct(read, d, p->getOverlapType(read, ud,p));
				for (int v=0; v<re.getVariationCount(d); v++)
					if (re.isMismatch(d, v))
						reporter.reportMismatch(v, getOverlapType(read.getData().isPositionInOverlap(d, read.getData().getMismatchPos(d, v)),read.getData(),d),retained.getQuick(v));
			}
					
			
		}
		fac.makeDistinct();
		SubreadsAlignedReadsData data = fac.createSubread();
		
		if (usedTotal!=null) usedTotal.accept(1, 1);
		ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> sread = new ImmutableReferenceGenomicRegion<>(sense?read.getReference():read.getReference().toOppositeStrand(), read.getRegion(), data);
		if (debug) {
			System.out.println("Convert:");
			System.out.println(read);
			System.out.println(sread);
			System.out.println();
		}
		
		return sread;
	}


	public static boolean getOverlapType(ImmutableReferenceGenomicRegion<? extends DefaultAlignedReadsData> read, int d,
			int p) {
		boolean in1 = !read.getData().hasGeometry() || read.getData().isPositionInFirstRead(d, p);
		boolean in2 = read.getData().hasGeometry() && read.getData().isPositionInSecondRead(d, p);
		return getOverlapType(in1&&in2, read.getData(), d);
	}

	public static boolean getOverlapType(boolean inBoth, AlignedReadsData data, int distinct) {
		if (data.getGeometryOverlap(distinct)==0) return false;
		if (!inBoth) return false;
		if (data.getGeometryBeforeOverlap(distinct)==0 && data.getGeometryAfterOverlap(distinct)==0) return true;
		return true;
	}

	
	public static void main(String[] args) {
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
		fac.start();
		fac.newDistinctSequence();
		fac.setGeometry(20, 100, 10);
		fac.addMismatch(10, 'C', 'T', false);
		ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read = ImmutableReferenceGenomicRegion.parse("1+:0-130",fac.create());
		System.out.println(read);
		System.out.println(new PairedEndToSubreadsConverter(false).convert(read, true, null,null));
		System.out.println(new PairedEndToSubreadsConverter(false).convert(read, false, null,null));
		
	}
	
}
