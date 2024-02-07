package gedi.grand10x;


import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public enum ReadCategory {
	
	
	Spliced,Unclear,Flank,FlankAntisense,Unspliced,Antisense,Intergenic;
	
	public static final int FLANK = 10_000;

	public static ReadCategory classify(Genomic g, ReferenceSequence ref, TrimmedGenomicRegion region) {
		boolean unspliced = false;
		boolean antisense = false;
		boolean flank = false;
		boolean flankantisense = false;
		
		for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(ref, region).loop()) {
			int exonlen = t.getRegion().intersect(region).getTotalLength();
			int intronlen = t.getRegion().invert().intersect(region).getTotalLength();
			if (exonlen>intronlen*10)
				return ReadCategory.Spliced;
			unspliced = true;
		}
		
		if (!unspliced)
			for (ImmutableReferenceGenomicRegion<String> tr : g.getGenes().ei(ref, region).list()) {
				unspliced = true;
			}
		
		for (ImmutableReferenceGenomicRegion<String> tr : g.getGenes().ei(ref.toOppositeStrand(), region).list()) {
			antisense = true;
		}

		if (!unspliced)
			for (ImmutableReferenceGenomicRegion<String> tr : g.getGenesWithFlank(FLANK).ei(ref, region).list()) {
				flank = true;
			}

		if (!antisense)
			for (ImmutableReferenceGenomicRegion<String> tr : g.getGenesWithFlank(FLANK).ei(ref.toOppositeStrand(), region).list()) {
				flankantisense = true;
			}

		
		for (ImmutableReferenceGenomicRegion<String> tr : g.getGenes().ei(ref.toOppositeStrand(), region).list()) {
			antisense = true;
		}

		if (unspliced && (antisense || flank || flankantisense)) return Unclear;
		if (unspliced) return Unspliced;
		
		if (antisense && (flank || flankantisense)) return Unclear;
		if (antisense) return Antisense;
		
		if (flank && flankantisense) return Unclear;
		if (flank) return Flank;
		
		if (flankantisense) return FlankAntisense;
		
		return ReadCategory.Intergenic;
		
	}
	
}
