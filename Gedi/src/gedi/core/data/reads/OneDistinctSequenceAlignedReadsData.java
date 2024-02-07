package gedi.core.data.reads;

import gedi.core.data.HasConditions;
import gedi.util.io.randomaccess.BinaryReader;

import java.io.IOException;

@Deprecated
public class OneDistinctSequenceAlignedReadsData extends SelectDistinctSequenceAlignedReadsData {

	public OneDistinctSequenceAlignedReadsData(AlignedReadsData parent,
			int distinct) {
		super(parent,distinct);
	}
		
}
