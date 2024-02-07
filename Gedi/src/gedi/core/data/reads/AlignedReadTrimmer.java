package gedi.core.data.reads;

import java.util.Iterator;
import java.util.function.Function;

import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public interface AlignedReadTrimmer<A extends AlignedReadsData> extends Function<ImmutableReferenceGenomicRegion<A>, Iterator<ImmutableReferenceGenomicRegion<A>>>{

	public AlignedReadTrimmer<A> setDebug(boolean debug);
	
	public static AlignedReadsVariation transformVariation(AlignedReadsVariation var, int start, int end, GenomicRegion mapper) {
		if (var.isSoftclip()) {
			return null;
		}
		
		// a deletion might be at the edge, and must then be trimmed!
		if(var.isDeletion()) {
			int delStart = var.getPosition();
			int delEnd = var.getPosition()+var.getReferenceSequence().length();
			if (delEnd<=start || delStart>=end)
				return null;
			
			int trimStart = Math.max(0, start-delStart);
			int trimEnd = Math.max(0, delEnd-end);
			
			int npos = mapper.induce(delStart+trimStart);
			CharSequence dels = var.getReferenceSequence().subSequence(trimStart, var.getReferenceSequence().length()-trimEnd);
			
			if (npos!=var.getPosition() || trimStart+trimEnd>0)
				return new AlignedReadsDeletion(npos, dels, var.isFromSecondRead());
			return var;
		}
		
		if (var.getPosition()<start || var.getPosition()>=end)
			return null;
		
		int npos = mapper.induce(var.getPosition());
		
		if (var.isMismatch()) {
			if (npos!=var.getPosition())
				return new AlignedReadsMismatch(npos, var.getReferenceSequence(), var.getReadSequence(), var.isFromSecondRead());
			return var;
		}
		else if(var.isInsertion()) {
			if (npos!=var.getPosition())
				return new AlignedReadsInsertion(npos, var.getReadSequence(), var.isFromSecondRead());
			return var;
		}
		throw new RuntimeException("Variant not supported!");
	}
}
