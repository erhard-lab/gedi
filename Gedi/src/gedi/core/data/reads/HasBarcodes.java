package gedi.core.data.reads;

import gedi.util.sequence.DnaSequence;

public interface HasBarcodes {

	
	 DnaSequence[] getBarcodes(int distinct, int condition);
	 int getBarcodeLength();
	 
}
