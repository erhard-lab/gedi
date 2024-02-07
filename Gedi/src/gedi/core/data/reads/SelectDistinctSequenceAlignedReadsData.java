package gedi.core.data.reads;

import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;

public class SelectDistinctSequenceAlignedReadsData implements AlignedReadsDataDecorator {

	private AlignedReadsData parent;
	private int[] distincts;
	
	public SelectDistinctSequenceAlignedReadsData(AlignedReadsData parent,
			int... distincts) {
		this.parent = parent;
		this.distincts = distincts;
	}
	
	public void setParent(AlignedReadsData parent) {
		this.parent = parent;
		hash = -1;
	}
	
	public AlignedReadsData getParent() {
		return parent;
	}
	
	@Override
	public boolean isFalseIntron(int l, int d) {
		return parent.isFalseIntron(l, this.distincts[d]);
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
		return parent.getGeometryBeforeOverlap(distincts[distinct]);
	}

	@Override
	public int getGeometryOverlap(int distinct) {
		return parent.getGeometryOverlap(distincts[distinct]);
	}

	@Override
	public int getGeometryAfterOverlap(int distinct) {
		return parent.getGeometryAfterOverlap(this.distincts[distinct]);
	}
	@Override
	public int getRawGeometry(int distinct) {
		return parent.getRawGeometry(distinct);
	}

	@Override
	public float getWeight(int distinct) {
		return parent.getWeight(this.distincts[distinct]);
	}
	
	@Override
	public int getId(int distinct) {
		return parent.getId(this.distincts[distinct]);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		throw new UnsupportedOperationException();
	}
	
	public void setDistinct(int distinct) {
		if (this.distincts.length==1)
			this.distincts[0] = distinct;
		else
			this.distincts = new int[] {distinct};
		hash = -1;
	}
	public void setDistinct(int... distinct) {
		this.distincts = distinct;
		hash = -1;
	}

	@Override
	public int getDistinctSequences() {
		return this.distincts.length;
	}

	@Override
	public int getNumConditions() {
		return parent.getNumConditions();
	}

	@Override
	public int getCount(int distinct, int condition) {
		return parent.getCount(this.distincts[distinct], condition);
	}

	@Override
	public int getVariationCount(int distinct) {
		return parent.getVariationCount(this.distincts[distinct]);
	}

	@Override
	public boolean isVariationFromSecondRead(int distinct, int index) {
		return parent.isVariationFromSecondRead(this.distincts[distinct], index);
	}

	
	@Override
	public boolean isMismatch(int distinct, int index) {
		return parent.isMismatch(this.distincts[distinct], index);
	}

	@Override
	public int getMismatchPos(int distinct, int index) {
		
		return parent.getMismatchPos(this.distincts[distinct], index);
	}

	@Override
	public CharSequence getMismatchGenomic(int distinct, int index) {
		
		return parent.getMismatchGenomic(this.distincts[distinct], index);
	}

	@Override
	public CharSequence getMismatchRead(int distinct, int index) {
		
		return parent.getMismatchRead(this.distincts[distinct], index);
	}

	@Override
	public boolean isInsertion(int distinct, int index) {
		
		return parent.isInsertion(this.distincts[distinct], index);
	}

	@Override
	public int getInsertionPos(int distinct, int index) {
		
		return parent.getInsertionPos(this.distincts[distinct], index);
	}

	@Override
	public CharSequence getInsertion(int distinct, int index) {
		
		return parent.getInsertion(this.distincts[distinct], index);
	}

	@Override
	public boolean isDeletion(int distinct, int index) {
		
		return parent.isDeletion(this.distincts[distinct], index);
	}

	@Override
	public int getDeletionPos(int distinct, int index) {
		
		return parent.getDeletionPos(this.distincts[distinct], index);
	}

	@Override
	public CharSequence getDeletion(int distinct, int index) {
		
		return parent.getDeletion(this.distincts[distinct], index);
	}
	
	@Override
	public boolean isSoftclip(int distinct, int index) {
		
		return parent.isSoftclip(this.distincts[distinct], index);
	}

	@Override
	public CharSequence getSoftclip(int distinct, int index) {
		
		return parent.getSoftclip(this.distincts[distinct], index);
	}

	@Override
	public boolean isSoftclip5p(int distinct, int index) {
		
		return parent.isSoftclip5p(this.distincts[distinct], index);
	}	
	

	@Override
	public int getMultiplicity(int distinct) {
		
		return parent.getMultiplicity(this.distincts[distinct]);
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
