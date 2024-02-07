package gedi.core.region;

import gedi.core.data.annotation.Transcript;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;

public enum GenomicRegionPosition {
	Start{
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			if (offset<0) return region.getStart()+offset;
			if (offset>region.getTotalLength()) return region.getEnd()+offset-region.getTotalLength();
			return region.mapMaybeOutside(offset);
		}
	},
	Stop {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			return End.position(reference, region, offset-1);
		}
	},
	End {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			if (offset>=0) return region.getEnd()+offset;
			if (offset<-region.getTotalLength()) return region.getStart()+offset+region.getTotalLength();
			return region.mapMaybeOutside(region.getTotalLength()+offset);
		}
	},
	FivePrime {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			return reference.getStrand()==Strand.Minus?Stop.position(reference,region,-offset):Start.position(reference,region,offset);
		}
	},
	ThreePrime {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			return reference.getStrand()==Strand.Minus?Start.position(reference,region,-offset):Stop.position(reference,region,offset);
		}
	}, 
	/** 
	 * If length is odd, the center is unique (e.g. 5 bp, the center is at index 2). If it is even, the center
	 * is the base closer to the 5' end. If region is spliced, the mapped induced position is taken!
	 */
	Center {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			int s = region.getTotalLength();
			if ((s&1)==1) s = s/2;
			else s=reference.getStrand()==Strand.Minus?s/2-1:s/2;
			return region.map(s)+(reference.getStrand()==Strand.Minus?offset:-offset);
		}
	},
	
	StartCodon {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			throw new RuntimeException("Must be called with referencesequence!");
		}
		
		@Override
		public int position(ReferenceGenomicRegion<?> ref, int offset) {
			if (!(ref.getData() instanceof Transcript)) throw new RuntimeException("Can only be called with Transcript data!");
			Transcript t = (Transcript) ref.getData();
			if (!t.isCoding()) throw new RuntimeException("Not a coding trancript!");
			return FivePrime.position(t.getCds(ref),offset);
		}
		
		@Override
		public boolean isValidInput(ReferenceGenomicRegion<?> ref) {
			if (!(ref.getData() instanceof Transcript))return false;
			Transcript t = (Transcript) ref.getData();
			if (!t.isCoding()) return false;
			return true;
		}
	},
	StopCodon {
		@Override
		public int position(ReferenceSequence reference, GenomicRegion region, int offset) {
			throw new RuntimeException("Must be called with referencesequence!");
		}
		
		@Override
		public int position(ReferenceGenomicRegion<?> ref, int offset) {
			if (!(ref.getData() instanceof Transcript)) throw new RuntimeException("Can only be called with Transcript data!");
			Transcript t = (Transcript) ref.getData();
			if (!t.isCoding()) throw new RuntimeException("Not a coding trancript!");
			return ThreePrime.position(t.getCds(ref),offset-2);
		}
		
		@Override
		public boolean isValidInput(ReferenceGenomicRegion<?> ref) {
			if (!(ref.getData() instanceof Transcript))return false;
			Transcript t = (Transcript) ref.getData();
			if (!t.isCoding()) return false;
			return true;
		}
	}
	
	;
	
	/**
	 * Offset is added depending on the value (start, stop w.r.t. genomic position, to 5' to 3' direction otherwise)
	 * @param reference
	 * @param region
	 * @param offset
	 * @return
	 */
	public abstract int position(ReferenceSequence reference, GenomicRegion region, int offset);
	
	public int position(ReferenceGenomicRegion<?> ref, int offset) {
		return position(ref.getReference(), ref.getRegion(), offset);
	}
	
	public boolean isValidInput(ReferenceGenomicRegion<?> ref) {
		return true;
	}
	
	public int position(ReferenceSequence reference, GenomicRegion region) {
		return position(reference, region, 0);
	}
	
	public int position(ReferenceGenomicRegion<?> ref) {
		return position(ref.getReference(), ref.getRegion(), 0);
	}
}