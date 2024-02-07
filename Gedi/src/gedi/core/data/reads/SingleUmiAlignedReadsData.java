package gedi.core.data.reads;

import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;

public class SingleUmiAlignedReadsData implements AlignedReadsData {

	private AlignedReadsData parent;
	private int distinct;
	private int condition;
	
	public SingleUmiAlignedReadsData(AlignedReadsData parent,
			int distinct, int condition) {
		this.parent = parent;
		this.distinct = distinct;
		this.condition = condition;
	}
	
	public int getOriginalCondition() {
		return condition;
	}
	
	public void setParent(AlignedReadsData parent) {
		this.parent = parent;
		hash = -1;
	}
	
	public AlignedReadsData getParent() {
		return parent;
	}
	
	@Override
	public boolean hasWeights() {
		return parent.hasWeights();
	}
	
	@Override
	public boolean hasGeometry() {
		return parent.hasGeometry();
	}

	@Override
	public int getGeometryBeforeOverlap(int distinct) {
		return parent.getGeometryBeforeOverlap(this.distinct);
	}

	@Override
	public int getGeometryOverlap(int distinct) {
		return parent.getGeometryOverlap(this.distinct);
	}

	@Override
	public int getGeometryAfterOverlap(int distinct) {
		return parent.getGeometryAfterOverlap(this.distinct);
	}
	
	@Override
	public float getWeight(int distinct) {
		return parent.getWeight(this.distinct);
	}
	
	@Override
	public int getRawGeometry(int distinct) {
		return parent.getRawGeometry(this.distinct);
	}

	@Override
	public int getId(int distinct) {
		return parent.getId(this.distinct);
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
	public int getCount(int distinct, int condition) {
		return 1;
	}
	
	@Override
	public boolean hasNonzeroInformation() {
		return true;
	}
	
	@Override
	public int getNonzeroCountValueForDistinct(int distinct, int index) {
		return 1;
	}
	@Override
	public int getNumConditionsWithCounts() {
		return 1;
	}
	
	@Override
	public int[] getNonzeroCountIndicesForDistinct(int distinct) {
		return new int[] {condition};
	}
	

	@Override
	public int getVariationCount(int distinct) {
		return parent.getVariationCount(this.distinct);
	}

	@Override
	public boolean isVariationFromSecondRead(int distinct, int index) {
		return parent.isVariationFromSecondRead(this.distinct, index);
	}

	
	@Override
	public boolean isMismatch(int distinct, int index) {
		return parent.isMismatch(this.distinct, index);
	}

	@Override
	public int getMismatchPos(int distinct, int index) {
		
		return parent.getMismatchPos(this.distinct, index);
	}

	@Override
	public CharSequence getMismatchGenomic(int distinct, int index) {
		
		return parent.getMismatchGenomic(this.distinct, index);
	}

	@Override
	public CharSequence getMismatchRead(int distinct, int index) {
		
		return parent.getMismatchRead(this.distinct, index);
	}

	@Override
	public boolean isInsertion(int distinct, int index) {
		
		return parent.isInsertion(this.distinct, index);
	}

	@Override
	public int getInsertionPos(int distinct, int index) {
		
		return parent.getInsertionPos(this.distinct, index);
	}

	@Override
	public CharSequence getInsertion(int distinct, int index) {
		
		return parent.getInsertion(this.distinct, index);
	}

	@Override
	public boolean isDeletion(int distinct, int index) {
		
		return parent.isDeletion(this.distinct, index);
	}

	@Override
	public int getDeletionPos(int distinct, int index) {
		
		return parent.getDeletionPos(this.distinct, index);
	}

	@Override
	public CharSequence getDeletion(int distinct, int index) {
		
		return parent.getDeletion(this.distinct, index);
	}
	
	@Override
	public boolean isSoftclip(int distinct, int index) {
		
		return parent.isSoftclip(this.distinct, index);
	}

	@Override
	public CharSequence getSoftclip(int distinct, int index) {
		
		return parent.getSoftclip(this.distinct, index);
	}

	@Override
	public boolean isSoftclip5p(int distinct, int index) {
		
		return parent.isSoftclip5p(this.distinct, index);
	}	
	

	@Override
	public int getMultiplicity(int distinct) {
		
		return parent.getMultiplicity(this.distinct);
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
