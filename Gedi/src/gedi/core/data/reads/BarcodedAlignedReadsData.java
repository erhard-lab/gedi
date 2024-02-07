package gedi.core.data.reads;

import java.io.IOException;
import java.util.Arrays;

import gedi.util.FileUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.sequence.DnaSequence;


public class BarcodedAlignedReadsData extends DefaultAlignedReadsData implements HasBarcodes {

	DnaSequence[][][] barcodes;
	

	public BarcodedAlignedReadsData() {}

	
	public DnaSequence[] getNonZeroBarcodes(int distinct, int index) {
		if (!hasNonzeroInformation()) throw new RuntimeException();
		return barcodes[distinct][index];
	}
	
	public DnaSequence[] getBarcodes(int distinct, int condition) {
		if (hasNonzeroInformation()) {
			
			int ind = Arrays.binarySearch(nonzeros[distinct], condition);
			if (ind<0) return new DnaSequence[0];
			return barcodes[distinct][ind];
		}
		return barcodes[distinct][condition];
	}
	
	
	public int getBarcodeLength() {
		for (int d=0; d<getDistinctSequences(); d++)
			for (int c=0; c<getNumConditions(); c++)
					if (barcodes[d][c].length>0) {
						return barcodes[d][c][0].length();
					}
		throw new RuntimeException("No barcode!");
	}
	
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		super.deserialize(in);
		int bl;
		DynamicObject gi = in.getContext().getGlobalInfo();
		if (!gi.hasProperty(BARCODEATTRIBUTE))
			bl = in.getCInt();
		else
			bl = gi.getEntry(BARCODEATTRIBUTE).asInt();
		
		if (hasNonzeroInformation()) {
			barcodes = new DnaSequence[getDistinctSequences()][][];
			
			for (int d=0; d<getDistinctSequences(); d++){
				barcodes[d] = new DnaSequence[count[d].length][];
				
				for (int i=0; i<count[d].length; i++) {
					barcodes[d][i] = new DnaSequence[count[d][i]];
					for (int b=0; b<count[d][i]; b++) {
						barcodes[d][i][b] = new DnaSequence(bl);
						FileUtils.readBitVector(barcodes[d][i][b], in,bl*2);
					}
				}
				
			}
			
		}
		else {
			barcodes = new DnaSequence[getDistinctSequences()][getNumConditions()][];
			for (int d=0; d<getDistinctSequences(); d++)
				for (int c=0; c<getNumConditions(); c++) {
					barcodes[d][c] = new DnaSequence[getCount(d, c)];
					for (int b=0; b<barcodes[d][c].length;b++) {
						barcodes[d][c][b] = new DnaSequence(bl);
						FileUtils.readBitVector(barcodes[d][c][b], in,bl*2);
					}
				}
		}
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (int d=0; d<getDistinctSequences(); d++) {
			if (hasId()) 
				sb.append(getId(d)).append(": ");
			if (hasNonzeroInformation()) {
				sb.append("[");
				int[] inds = getNonzeroCountIndicesForDistinct(d);
				for (int i=0; i<inds.length; i++)
					sb.append(inds[i]+":"+getNonzeroCountValueForDistinct(d, i)+",");
				sb.deleteCharAt(sb.length()-1);
				sb.append("]");
			} else
				sb.append(Arrays.toString(getCountsForDistinctInt(d, ReadCountMode.All)));
			sb.append(" x");
			sb.append(getMultiplicity(d));
			if (hasWeights()) 
				sb.append(" (w=").append(String.format("%.2f", getWeight(d))).append(")");
			if (hasGeometry()) 
				sb.append(String.format(" %d|%d|%d", getGeometryBeforeOverlap(d),getGeometryOverlap(d),getGeometryAfterOverlap(d)));
			for (AlignedReadsVariation var : getVariations(d))
				sb.append("\t"+var);
			sb.append(" [");
			for (int c=0; c<getNumConditions(); c++) {
				if (c>0) sb.append(",");
				DnaSequence[] bcs = getBarcodes(d, c);
				for (int b=0; b<bcs.length; b++) {
					if (b>0) sb.append(" ");
					sb.append(bcs[b]);
				}
			}
			sb.append("]");
			if (d<getDistinctSequences()-1) sb.append(" ~ ");
		}
		return sb.toString();
	}



	

}
