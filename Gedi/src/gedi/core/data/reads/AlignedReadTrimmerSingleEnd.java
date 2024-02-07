package gedi.core.data.reads;

import java.util.TreeMap;

import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

public class AlignedReadTrimmerSingleEnd<A extends AlignedReadsData> implements AlignedReadTrimmer<A> {

	private int start;
	private int end;
	
	public AlignedReadTrimmerSingleEnd(int start, int end) {
		this.start = start;
		this.end = end;
	}
	
	private boolean debug = false;
	
	public AlignedReadTrimmerSingleEnd<A> setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	
	@Override
	public ExtendedIterator<ImmutableReferenceGenomicRegion<A>> apply(
			ImmutableReferenceGenomicRegion<A> t) {
		
		if (t.getData().hasGeometry()) throw new RuntimeException("Cannot clip a paired end read with AlignedReadClipperSingleEnd!");
		
		int len = t.getRegion().getTotalLength();
		
		TreeMap<GenomicRegion,AlignedReadsDataFactory> map = new TreeMap<GenomicRegion, AlignedReadsDataFactory>();
		
		for (int d=0; d<t.getData().getDistinctSequences(); d++) {
			
			GenomicRegion r = new ArrayGenomicRegion(0,len);
			
			// adapt the trimming if read is softclipped
			int start = this.start;
			int end = this.end;
			for (int v=0; v<t.getData().getVariationCount(d); v++) 
				if (t.getData().isSoftclip(d, v)) {
					int l = t.getData().getSoftclip(d, v).length();
					if (t.getData().isSoftclip5p(d, v))
						start-=l;
					else
						end-=l;
				}
			
			// adapt the trimming if the desired trimming is within a deletion
			for (int v=0; v<t.getData().getVariationCount(d); v++) 
				if (t.getData().isDeletion(d, v)) {
					int dstart = t.getData().getDeletionPos(d, v);
					int dend = dstart+t.getData().getDeletion(d, v).length();
					if (start>=dstart && start<dend) 
						start = dend;
					else if (r.getTotalLength()-end>dstart && r.getTotalLength()-end<=dend)
						end=r.getTotalLength()-dstart;
				}
			start = Math.max(start, 0);
			end = Math.max(end, 0);
			
			GenomicRegion treg = r.map(start, r.getTotalLength()-end); // is empty if length would be negative!
			// this is still in the induced coords system of t!
			if (treg.isEmpty()) 
				continue; // everything is clipped, might happen for very short fragments!
			AlignedReadsDataFactory fac = map.computeIfAbsent(t.map(treg), x-> 
						new AlignedReadsDataFactory(t.getData().getNumConditions(),t.getData().hasNonzeroInformation()).start());
			
			fac.add(t.getData(), d, v->AlignedReadTrimmer.transformVariation(v,treg.getStart(),treg.getEnd(),treg),false);
			
			
		}
		
		Class<A> cls = (Class<A>) t.getData().getClass();
		if (debug) {
			System.out.println("Trim:");
			System.out.println(t);
			EI.wrap(map.keySet()).map(reg->{
				AlignedReadsDataFactory fac = map.get(reg);
				fac.makeDistinct();
				return new ImmutableReferenceGenomicRegion<>(t.getReference(), reg, fac.create(cls));
			}).print();
			return EI.wrap(map.keySet()).map(reg->{
				AlignedReadsDataFactory fac = map.get(reg);
				fac.makeDistinct();
				return new ImmutableReferenceGenomicRegion<>(t.getReference(), reg, fac.create(cls));
			});
		}
		
		return EI.wrap(map.keySet()).map(reg->{
			AlignedReadsDataFactory fac = map.get(reg);
			fac.makeDistinct();
			return new ImmutableReferenceGenomicRegion<>(t.getReference(), reg, fac.create(cls));
		});

	}
	

}
