package gedi.grand10x;


import java.util.ArrayList;

import gedi.core.data.annotation.Transcript;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;

/**
 * Important for Cell Ranger mapped reads, as these might go into an intron by 3 nt!
 * @author erhard
 *
 */
public class TrimmedGenomicRegion implements GenomicRegion {
	private GenomicRegion parent;
	
	public TrimmedGenomicRegion set(GenomicRegion reg) {
		this.parent = reg;
		return this;
	}
	
	public GenomicRegion getParent() {
		return parent;
	}
	
	@Override
	public String toString() {
		return toString2();
	}
	
	@Override
	public int getStart(int part) {
		if (!trimStart()) 
			return part==0?parent.getStart()+3:parent.getStart(part);
		return parent.getStart(part+1);
	}
	
	@Override
	public int getEnd(int part) {
		if (part==getNumParts()-1)
			return trimEnd()?parent.getEnd(parent.getNumParts()-2):parent.getEnd()-3;
			
		return parent.getEnd(trimStart()?part+1:part);
	}
	
	private final boolean trimStart() {
		return parent.getLength(0)<=3;
	}
	
	private final boolean trimEnd() {
		return parent.getLength(parent.getNumParts()-1)<=3;
	}
	
	@Override
	public int getNumParts() {
		return parent.getNumParts()-(trimStart()?1:0)-(trimEnd()?1:0);
	}
	
	
	public boolean isCompatibleWith(GenomicRegion t) {
		return t.intersect(this).getNumParts()==this.getNumParts() 
				&& t.isIntronConsistent(this);
	}

	
}