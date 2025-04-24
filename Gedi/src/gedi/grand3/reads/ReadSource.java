package gedi.grand3.reads;

import gedi.core.data.reads.AlignedReadTrimmer;
import gedi.core.data.reads.AlignedReadTrimmerPairedEnd;
import gedi.core.data.reads.AlignedReadTrimmerSingleEnd;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.data.reads.subreads.MismatchReporter;
import gedi.core.data.reads.subreads.ToSubreadsConverter;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.ExtendedIterator;


/**
 * Iterate over all reads from a location while do on-the-fly trimming. What "from a location" means depends on the given strandness 
 * @author erhard
 *
 * @param <A>
 */
public class ReadSource<A extends AlignedReadsData> {
	
	private AlignedReadTrimmer<A> trimmer;
	private Strandness strandness;
	public GenomicRegionStorage<A> reads;
	private ToSubreadsConverter<A> converter;
	
	public ReadSource(GenomicRegionStorage<A> reads, ClippingData clipping, Strandness strandness, boolean debug) {
		this.reads = reads;
		this.strandness = strandness;
		
		if (clipping!=null && !clipping.isNoClipping()) {
			if (reads.getRandomRecord().hasGeometry())
				trimmer = new AlignedReadTrimmerPairedEnd<A>(clipping.getClip5p1(), clipping.getClip5p2(), clipping.getClip3p1(), clipping.getClip3p2()).setDebug(debug);
			else
				trimmer = new AlignedReadTrimmerSingleEnd<A>(clipping.getClip5p1(),clipping.getClip3p1()).setDebug(debug);
			trimmer.setDebug(debug);
		} else
			trimmer = null;
		
		converter = ToSubreadsConverter.infer(reads, debug);
	}
	
	public ToSubreadsConverter<A> getConverter() {
		return converter;
	}
	
	public Strandness getStrandness() {
		return strandness;
	}
	
	public ExtendedIterator<ImmutableReferenceGenomicRegion<A>> getRawReads(ImmutableReferenceGenomicRegion<?> location){

		if (strandness.equals(Strandness.Sense))
			return reads.ei(location);
		if (strandness.equals(Strandness.Antisense))
			return reads.ei(location.toOppositeStrand());

			return reads.ei(location).chain(reads.ei(location.toOppositeStrand()));
	}
	

	public ExtendedIterator<ImmutableReferenceGenomicRegion<A>> getReads(ImmutableReferenceGenomicRegion<?> location){

		if (strandness.equals(Strandness.Sense))
			return reads.ei(location).iff(trimmer!=null, ei->ei.unfold(trimmer));
		if (strandness.equals(Strandness.Antisense))
			return reads.ei(location.toOppositeStrand()).iff(trimmer!=null, ei->ei.unfold(trimmer));

			return reads.ei(location).chain(reads.ei(location.toOppositeStrand())).iff(trimmer!=null, ei->ei.unfold(trimmer));
	}
	
	
	public ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> getSubReads(ImmutableReferenceGenomicRegion<String> location, MismatchReporter reporter){
		return converter.convert(location.toLocationString()+":"+location.getData(),location.getReference(), getReads(location).list(), reporter,null);
	}
	
}
