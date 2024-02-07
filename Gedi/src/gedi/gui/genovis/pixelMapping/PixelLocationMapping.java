package gedi.gui.genovis.pixelMapping;

import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.util.gui.PixelBasepairMapper;

import java.util.ArrayList;

/**
 * Guarranteed that blocks are sorted according to their screen coordinates from left to right
 * @author erhard
 *
 */
public class PixelLocationMapping extends ArrayList<PixelLocationMappingBlock> {

	
	private PixelBasepairMapper xmapper;
	
	
	public PixelLocationMapping() {
		
	}
	
	public PixelLocationMapping(ReferenceSequence reference, PixelBasepairMapper xmapper) {
		this.xmapper = xmapper;
		GenomicRegion region = xmapper.getRegion(reference);
		for (int p=0; p<region.getNumParts(); p++) {
			for (int i=region.getStart(p); i<region.getEnd(p); i++) {
				int start = i;
				double s = xmapper.bpToPixel(reference,i);
				double e = xmapper.bpToPixel(reference,i+1);
				if (Math.abs(e-s)<1) {
					int x = xmapper.pixelToBp(s+1, reference);
					i = x==-1 || region.getStop(p)<x?region.getStop(p):x;
					e = xmapper.bpToPixel(reference,i+1);
				}
				PixelLocationMappingBlock mb = obtainMappingBlock();
				mb.set(reference, start, i);
			}
		}
		
		if (xmapper.is5to3() && reference.getStrand()==Strand.Minus) {
			int from = 0;
			int to = size()-1;
	        while (from< to) {
	        	set(from++,set(to++,get(from)));
	        }
		}
	}
	
	public int getBlockForBp(ReferenceSequence reference, int bp) {
		if (size()==0) return -1;
		double pix = xmapper.bpToPixel(reference, bp);
		double rel = pix/xmapper.getWidth();
		int re = (int) (rel*size());
//		if (re==0 && !get(re).containsBp(bp))
//			return -1;
		if (re<size() && re>=0 && !get(re).containsBp(bp)) {
			if (re>0 && get(re-1).containsBp(bp)) return re-1;
			if (re+1<size() && get(re+1).containsBp(bp)) return re+1;
			return binarySearch(bp);
		}
		return re;
	}
	
		
	private int binarySearch(int bp) {
		int low = 0;
        int high = size()-1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            PixelLocationMappingBlock midVal = get(mid);
            int cmp = midVal.compareToBp(bp);

            if (cmp < 0)
                low = mid + 1;
            else if (cmp > 0)
                high = mid - 1;
            else
                return mid; // key found
        }
        if (high==-1) return -1;
        if (low==size()) return size();
        throw new RuntimeException("nononono");
	}

	public PixelBasepairMapper getXmapper() {
		return xmapper;
	}

	private PixelLocationMappingBlock obtainMappingBlock() {
		PixelLocationMappingBlock re = new PixelLocationMappingBlock();
		add(re);
		return re;
	}
	
	public PixelLocationMappingBlock addBlock(ReferenceSequence ref, int startBp, int stopBp) {
		PixelLocationMappingBlock bl = obtainMappingBlock();
		bl.set(ref, startBp, stopBp);
		return bl;
	}
	
	@Override
	public String toString() {
		return getRegion().toString();
	}



	public GenomicRegion getRegion() {
		ArrayGenomicRegion re = new ArrayGenomicRegion();
		for (PixelLocationMappingBlock b : this)
			re = re.union(new ArrayGenomicRegion(b.getStartBp(),b.getStopBp()+1));
		return re;
	}
	
}

