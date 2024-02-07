package gedi.core.sequence;

import java.util.Set;

import gedi.core.data.annotation.ReferenceSequenceLengthProvider;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;

public interface SequenceProvider extends ReferenceSequenceLengthProvider {

	
	/**
	 * Instead of throwing exceptions, returns the appropriate amount of N, if sequence not available (i.e. ref not known, region exceeds bounds).
	 * @param rgr
	 * @return
	 */
	default CharSequence getSequenceSave(ReferenceGenomicRegion<?> rgr) {
		if (rgr==null) return null;
		if (!knowsSequence(rgr.getReference().getName()) || !rgr.getRegion().intersects(getRegionOfReference(rgr.getReference().getName()))) return StringUtils.repeatSequence('N', rgr.getRegion().getTotalLength());

		CharSequence re = getSequence(rgr.getReference(),rgr.getRegion().intersect(getRegionOfReference(rgr.getReference().getName())));
		
		if (rgr.getRegion().getStart()<0) re = StringUtils.repeatSequence('N', new ArrayGenomicRegion(rgr.getRegion().getStart(),0).intersect(rgr.getRegion()).getTotalLength())+re.toString();
		if (rgr.getRegion().getEnd()>getLength(rgr.getReference().getName())) re = re.toString()+StringUtils.repeatSequence('N', new ArrayGenomicRegion(getLength(rgr.getReference().getName()),rgr.getRegion().getEnd()).intersect(rgr.getRegion()).getTotalLength());
		
		return re;
	}
	
	default CharSequence getSequenceSave(ReferenceSequence ref, GenomicRegion region) {
		if (ref==null || region==null) return null;
		if (!knowsSequence(ref.getName()) || !region.intersects(getRegionOfReference(ref.getName()))) return StringUtils.repeatSequence('N', region.getTotalLength());

		CharSequence re = getSequence(ref,region.intersect(getRegionOfReference(ref.getName())));
		
		if (region.getStart()<0) re = StringUtils.repeatSequence('N', new ArrayGenomicRegion(region.getStart(),0).intersect(region).getTotalLength())+re.toString();
		if (region.getEnd()>getLength(ref.getName())) re = re.toString()+StringUtils.repeatSequence('N', new ArrayGenomicRegion(getLength(ref.getName()),region.getEnd()).intersect(region).getTotalLength());
		
		return re;
	}
	
	
	/**
	 * In 5' to 3' direction! Returns null if ref is unknown; throws an Exception if region is outside of ref's bounds
	 * @param ref
	 * @param region
	 * @return
	 */
	default CharSequence getSequence(ReferenceGenomicRegion<?> rgr) {
		if (rgr==null) return null;
		return getSequence(rgr.getReference(),rgr.getRegion());
	}
	
	
	/**
	 * In 5' to 3' direction! Returns null if ref is unknown; throws an Exception if region is outside of ref's bounds
	 * @param ref
	 * @param region
	 * @return
	 */
	default CharSequence getSequence(ReferenceSequence ref, GenomicRegion region) {
		CharSequence re = getPlusSequence(ref.getName(),region);
		if (re==null) return null;
		if (ref.getStrand()==Strand.Minus)
			re = SequenceUtils.getDnaReverseComplement(re);
		return re;
	}
	default char getSequence(ReferenceSequence ref, int pos) {
		char re = getPlusSequence(ref.getName(),pos);
		if (ref.getStrand()==Strand.Minus)
			re = SequenceUtils.getDnaComplement(re);
		return re;
	}
	default boolean knowsSequence(String name) {
		return getSequenceNames().contains(name);
	}
	
	
	Set<String> getSequenceNames();
	CharSequence getPlusSequence(String name, GenomicRegion region);
	char getPlusSequence(String name, int pos);
	
	
}
