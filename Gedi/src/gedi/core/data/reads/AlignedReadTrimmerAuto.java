package gedi.core.data.reads;

import java.io.IOException;

import gedi.app.Gedi;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.ExtendedIterator;

public class AlignedReadTrimmerAuto<A extends AlignedReadsData> implements AlignedReadTrimmer<A> {

	private AlignedReadTrimmerSingleEnd<A> single;
	private AlignedReadTrimmerPairedEnd<A> paired;
	
	public AlignedReadTrimmerAuto(int start1, int start2, int end1, int end2) {
		this.single = new AlignedReadTrimmerSingleEnd<>(start1, end1);
		this.paired = new AlignedReadTrimmerPairedEnd<>(start1, start2, end1, end2);
	}
	
	
	public AlignedReadTrimmerAuto<A> setDebug(boolean debug) {
		single.setDebug(debug);
		paired.setDebug(debug);
		return this;
	}
	
	
	@Override
	public ExtendedIterator<ImmutableReferenceGenomicRegion<A>> apply(
			ImmutableReferenceGenomicRegion<A> t) {
		
		if (t.getData().hasGeometry()) 
			return paired.apply(t);
		else
			return single.apply(t);

	}
	
	
	public static void main(String[] args) throws IOException {
		Gedi.startup(false);
		GenomicRegionStorage<AlignedReadsData> st = Gedi.load(args[0]);
		AlignedReadTrimmerAuto<AlignedReadsData> c = new AlignedReadTrimmerAuto<>(25,25,35,15);

		for (ImmutableReferenceGenomicRegion<AlignedReadsData> r : st.ei().loop()) {
//			if (c.apply(r).count()!=1) {
				System.out.println(r);
				for (ImmutableReferenceGenomicRegion<AlignedReadsData> t : c.apply(r).loop()) 
					System.out.println(t);
				System.out.println();
//			}
		}
	}



}
