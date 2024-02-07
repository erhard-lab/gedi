package gedi.core.data.reads;

import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;

/**
 * You have to specify the read count mode in the constructor, such that this class knows what to do with several distinct sequences
 * and different weights!
 * 
 * If you use {@link ReadCountMode#Weight} to retrieve counts from this object, you will get the sum of what you would have gotten
 * for using the given mode on the parent's distinct sequences.
 * @author erhard
 *
 */
public class IgnoreVariationsAlignedReadsData implements AlignedReadsDataDecorator {

	private AlignedReadsData parent;
	private ReadCountMode mode;
	
	public IgnoreVariationsAlignedReadsData(AlignedReadsData parent, ReadCountMode mode) {
		this.parent = parent;
		this.mode = mode;
	}
	
	public void setParent(AlignedReadsData parent) {
		this.parent = parent;
		hash = -1;
	}
	
	public AlignedReadsData getParent() {
		return parent;
	}
	
	
	@Override
	public int getId(int distinct) {
		if (distinct!=0) throw new IndexOutOfBoundsException();
		return parent.getId(0);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		throw new UnsupportedOperationException();
	}
	

	@Override
	public int getDistinctSequences() {
		return 1;
	}

	@Override
	public int getNumConditions() {
		return parent.getNumConditions();
	}
	
	@Override
	public boolean hasGeometry() {
		return parent.hasGeometry();
	}

	@Override
	public int getGeometryBeforeOverlap(int distinct) {
		return parent.getGeometryBeforeOverlap(distinct);
	}

	@Override
	public int getGeometryOverlap(int distinct) {
		return parent.getGeometryOverlap(distinct);
	}

	@Override
	public int getGeometryAfterOverlap(int distinct) {
		return parent.getGeometryAfterOverlap(distinct);
	}
	
	@Override
	public int getRawGeometry(int distinct) {
		return parent.getRawGeometry(distinct);
	}

	@Override
	public boolean hasWeights() {
		return true;
	}
	
	@Override
	public float getWeight(int distinct) {
		if (distinct!=0) throw new IndexOutOfBoundsException();
		
		double mtot = 0;
		double alltot = 0;
		for (distinct=0; distinct<parent.getDistinctSequences(); distinct++) {
			mtot+=parent.getTotalCountForDistinct(distinct, mode);
			alltot+=parent.getTotalCountForDistinct(distinct, ReadCountMode.All);
		}
		
		return (float) (mtot/alltot);
	}
	
	@Override
	public int getCount(int distinct, int condition) {
		if (distinct!=0) throw new IndexOutOfBoundsException();
		return parent.getTotalCountForConditionInt(condition, ReadCountMode.All);
	}
	
	@Override
	public int getMultiplicity(int distinct) {
		return 0;
	}


	@Override
	public int getVariationCount(int distinct) {
		return 0;
	}

	@Override
	public boolean isVariationFromSecondRead(int distinct, int index) {
		return false;
	}

	
	@Override
	public boolean isMismatch(int distinct, int index) {
		return false;
	}

	@Override
	public int getMismatchPos(int distinct, int index) {
		return -1;
	}

	@Override
	public CharSequence getMismatchGenomic(int distinct, int index) {
		return null;
	}

	@Override
	public CharSequence getMismatchRead(int distinct, int index) {
		return null;
	}

	@Override
	public boolean isInsertion(int distinct, int index) {
		return false;
	}

	@Override
	public int getInsertionPos(int distinct, int index) {
		return -1;
	}

	@Override
	public CharSequence getInsertion(int distinct, int index) {
		return null;
	}

	@Override
	public boolean isDeletion(int distinct, int index) {
		return false;
	}

	@Override
	public int getDeletionPos(int distinct, int index) {
		return -1;
	}

	@Override
	public CharSequence getDeletion(int distinct, int index) {
		return null;
	}
	
	@Override
	public boolean isSoftclip(int distinct, int index) {
		return false;
	}

	@Override
	public CharSequence getSoftclip(int distinct, int index) {
		return null;
	}

	@Override
	public boolean isSoftclip5p(int distinct, int index) {
		return false;
	}	
	

	
	transient int hash = -1;
	@Override
	public int hashCode() {
		if (hash==-1) hash = hashCode2();
		return hash;
	}
	@Override
	public boolean equals(Object obj) {
		return equals(obj,true,true);
	}
	
	@Override
	public String toString() {
		return toString2();
	}


	
	
	
}
