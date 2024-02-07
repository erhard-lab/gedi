package gedi.core.data.reads.subreads;

import gedi.core.data.reads.AlignedReadsData;

public interface ReadByReadToSubreadsConverter<A extends AlignedReadsData> extends ToSubreadsConverter<A> {


	@Override
	default boolean isReadByRead() {
		return true;
	}
	
}