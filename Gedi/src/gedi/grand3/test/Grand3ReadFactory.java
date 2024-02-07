package gedi.grand3.test;

import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.functions.EI;
import gedi.util.math.stat.RandomNumbers;

public class Grand3ReadFactory {

	private Genomic genomic;
	private RandomNumbers rnd = new RandomNumbers(42);
	private int retries = 100;
	
	
	public Grand3ReadFactory(Genomic genomic) {
		this.genomic = genomic;
	}


	public ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> create(String loc, 
			boolean sense, boolean exon, 
			int geometryBefore, int geometryOverlap, int geometryAfter, 
			int kBefore, int kOverlap, int kInconsistentOverlap, int kAfter,
			int[] with, int[] without) {
		
		
		ImmutableReferenceGenomicRegion l = ImmutableReferenceGenomicRegion.parse(genomic, loc);
		l = genomic.getTranscripts().ei(l).filter(t->t.getRegion().getNumParts()>1).first();
		
		ReferenceSequence ref = sense?l.getReference():l.getReference().toOppositeStrand();
		
		int len = geometryBefore+geometryOverlap+geometryAfter;
		
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(with.length).start();
		fac.newDistinctSequence();
		fac.setCount(without);
		fac.setMultiplicity(1);
		fac.setGeometry(geometryBefore,geometryOverlap,geometryAfter);

		fac.newDistinctSequence();
		fac.setCount(with);
		fac.setMultiplicity(1);
		fac.setGeometry(geometryBefore,geometryOverlap,geometryAfter);
		
		int[] before = null;
		int[] overlap = null;
		int[] after = null;
		GenomicRegion region = null;
		
		for (int i=0; i<retries; i++) {
			if (exon) {
				int s = rnd.getUnif(0, l.getRegion().getTotalLength()-len+1);
				region = l.map(new ArrayGenomicRegion(s,s+len));
			}
			else {
				GenomicRegion intr = l.getRegion().invert();
				if (intr.isEmpty()) return null;
				intr = intr.getPart(rnd.getUnif(0, intr.getNumParts())).asRegion();
				int s = rnd.getUnif(0, intr.getTotalLength()+len-1);
				region = new ArrayGenomicRegion(intr.getStart()-len+1+s,intr.getStart()-len+1+s+len);
			}
	
			char[] seq = genomic.getSequence(ref,region).toString().toCharArray();
			before = getConvPos(seq, 0, geometryBefore, kBefore,sense?'T':'A');
			if (before==null) continue;
			overlap = getConvPos(seq, geometryBefore, geometryBefore+geometryOverlap, kOverlap+kInconsistentOverlap,sense?'T':'A');
			if (overlap==null) continue;
			after = getConvPos(seq, geometryBefore+geometryOverlap, geometryBefore+geometryOverlap+geometryAfter, kAfter,sense?'T':'A');
			if (after==null) continue;
			break;
		}
		
		if (before==null || overlap==null || after==null)
			return null;

		
		for (int pos : before)
			fac.addMismatch(pos, sense?'T':'A', sense?'C':'G', false);
		
		for (int pi=0; pi<kOverlap; pi++) {
			fac.addMismatch(overlap[pi], sense?'T':'A', sense?'C':'G', false);
			fac.addMismatch(overlap[pi], !sense?'T':'A', !sense?'C':'G', true);
		}
		for (int pi=kOverlap; pi<overlap.length; pi++) {
			fac.addMismatch(overlap[pi], sense?'T':'A', sense?'C':'G', false);
		}

		for (int pos : after)
			fac.addMismatch(pos, !sense?'T':'A', !sense?'C':'G', true);

		
		return new ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>(ref, region, fac.create());
	}
	
	
	private int[] getConvPos(char[] seq, int start, int end, int k, char g) {
		int[] tpos = EI.seq(start, end).filterInt(pos->seq[pos]==g).toIntArray();
		if (tpos.length<k) return null;
		rnd.shuffle(tpos);
		return ArrayUtils.slice(tpos, 0, k);
	}
	
}
