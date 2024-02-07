package gedi.core.data.reads.subreads;

import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public class NoopToSubreadsConverter implements ReadByReadToSubreadsConverter<SubreadsAlignedReadsData> {
	
	private String[] semantic;
	
	
	public NoopToSubreadsConverter(String[] semantic) {
		this.semantic = semantic;
	}

	@Override
	public String[] getSemantic() {
		return semantic;
	}
	
	@Override
	public ToSubreadsConverter setDebug(boolean debug) {
		return this;
	}
	
	
	@Override
	public ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> convert(
			ImmutableReferenceGenomicRegion<? extends SubreadsAlignedReadsData> read, boolean sense, MismatchReporter reporter) {
		
		SubreadsAlignedReadsData re = read.getData();
		for (int d=0; d<re.getDistinctSequences(); d++) {
			// report retained/corrected mismatches
			if (reporter!=null) {
				reporter.startDistinct(read, d, p->false);
				for (int v=0; v<re.getVariationCount(d); v++)
					if (re.isMismatch(d, v))
						reporter.reportMismatch(v, re.getSubreadIndexForPosition(d, re.getMismatchPos(d, v), read.getRegion().getTotalLength())>0, true);
			}
		}
		
		return (ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>) read;
		
	}

	
	

}
