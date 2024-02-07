package gedi.core.data.reads;

import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.datastructure.collections.intcollections.IntIterator;

public interface HasSubreads {
	
	 
	/**
	 * two reads might be joined, that do not overlap (e.g. for a UMI, or paired-end reads). In this case the region will have an "intron"
	 * that actually is no intron. To evaluate consistency with e.g. a transcript, this information is provided here. This returns
	 * the positions in the read where such "introns" start (###----### would be position 3)
	 * 
	 * Be careful: this end might not necessarily refer to a position, where an intron starts:
	 * ##########
	 *           ######---------####
	 * here, position 10 will be returned, which is the subread end 10.    
	 *          
	 * @return
	 */
	IntIterator getGapPositions(int distinct);
	int getNumSubreads(int distinct);
	
	int getSubreadStart(int distinct, int index);
	int getSubreadId(int distinct, int index);
	
	
	default boolean isConsistentlyContained(ReferenceGenomicRegion<?> read,
			ReferenceGenomicRegion<?> reference, int distinct) {
		int last = 0;
		IntIterator it = getGapPositions(distinct);
		while (it.hasNext()) {
			int c = it.nextInt();
			GenomicRegion part = read.map(new ArrayGenomicRegion(last,c));
			if (!reference.getRegion().containsUnspliced(part)) return false;
			last = c;
		}
		GenomicRegion part = last==0?read.getRegion():read.map(new ArrayGenomicRegion(last,read.getRegion().getTotalLength()));
		if (!reference.getRegion().containsUnspliced(part)) return false;
		
		return true;
	}
	
	default boolean isFalseIntron(int pos, int distinct) {
		return getGapPositions(distinct).filterInt(p->p==pos).count()>0;
	}
	default int getNumParts(ReferenceGenomicRegion<?> read,int distinct) {
		return read.getRegion().getNumParts()-
				getGapPositions(distinct)
				.filterInt(p->Math.abs(read.map(p)-read.map(p-1))!=1).countInt();
	}

	
	
	default int getSubreadIndexForPosition(int distinct, int position, int totalLength) {
		int low = 0;
        int high = getNumSubreads(distinct)-1;
        position++; // getSubreadEnd which is compared against is exclusive!
        while (low <= high) {
            int mid = (low + high) >>> 1;
            int cmp = Integer.compare(getSubreadEnd(distinct,mid,totalLength), position);

            if (cmp < 0)
                low = mid + 1;
            else if (cmp > 0)
                high = mid - 1;
            else
                return mid; // key found
        }
        return low;  // key not found.
	}
	
	
	
	default int getSubreadEnd(int distinct, int index, int totalLength) {
		if (index==getNumSubreads(distinct)-1) return totalLength;
		return getSubreadStart(distinct,index+1);
	}
	
	
	
	default StringBuilder addSubreadToString(StringBuilder sb, int d) {
		sb.append("<");
		for (int b=0; b<getNumSubreads(d); b++) {
			sb.append(getSubreadId(d,b));
			if (b<getNumSubreads(d)-1)
				sb.append("@").append(getSubreadEnd(d,b,0)).append("/");
		}
		sb.append(">");
		return sb;
	}
	
	
	
}
