package gedi.core.data.reads;

import java.util.TreeMap;

import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

public class AlignedReadTrimmerPairedEnd<A extends AlignedReadsData> implements AlignedReadTrimmer<A> {

	private int start1;
	private int start2;
	private int end1;
	private int end2;
	
	public AlignedReadTrimmerPairedEnd(int start1, int start2, int end1, int end2) {
		this.start1 = start1;
		this.start2 = start2;
		this.end1 = end1;
		this.end2 = end2;
	}
	
	private boolean debug = false;
	
	public AlignedReadTrimmerPairedEnd<A> setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	
	@Override
	public ExtendedIterator<ImmutableReferenceGenomicRegion<A>> apply(
			ImmutableReferenceGenomicRegion<A> t) {
		
		if (!t.getData().hasGeometry()) throw new RuntimeException("Cannot clip a single end read with AlignedReadClipperPairedEnd!");
		
		int len = t.getRegion().getTotalLength();
		
		TreeMap<GenomicRegion,AlignedReadsDataFactory> map = new TreeMap<GenomicRegion, AlignedReadsDataFactory>();

		for (int d=0; d<t.getData().getDistinctSequences(); d++) {
			
			int before = t.getData().getGeometryBeforeOverlap(d);
			int overlap = t.getData().getGeometryOverlap(d);
			int after = t.getData().getGeometryAfterOverlap(d);
			
			if (before+overlap+after!=len)
				throw new RuntimeException("Your mapping is broken, geometry does not match the covered region!");
			
			GenomicRegion r1 = new ArrayGenomicRegion(0,before+overlap);
			GenomicRegion r2 = new ArrayGenomicRegion(before,before+overlap+after);
			
			// adapt the trimming if read is softclipped
			int start1 = this.start1;
			int start2 = this.start2;
			int end1 = this.end1;
			int end2 = this.end2;
			for (int v=0; v<t.getData().getVariationCount(d); v++) 
				if (t.getData().isSoftclip(d, v)) {
					int l = t.getData().getSoftclip(d, v).length();
					if (t.getData().isVariationFromSecondRead(d, v) ) {
						if (t.getData().isSoftclip5p(d, v))
							end2-=l;
						else
							start2-=l;
					} else if (t.getData().isSoftclip5p(d, v))
						start1-=l;
					else
						end1-=l;
				}
			
			// adapt the trimming if the desired trimming is within a deletion
			for (int v=0; v<t.getData().getVariationCount(d); v++) 
				if (t.getData().isDeletion(d, v)) {
					int dstart = t.getData().getDeletionPos(d, v);
					int dend = dstart+t.getData().getDeletion(d, v).length();
					if (t.getData().isVariationFromSecondRead(d, v) ) {
						if (len-start2>dstart && len-start2<=dend) 
							start2 = len-dstart;
						else if (before+end2>=dstart && before+end2<dend)
							end2=dend-before;
					} else {
						if (start1>=dstart && start1<dend) 
							start1 = dend;
						else if (r1.getTotalLength()-end1>dstart && r1.getTotalLength()-end1<=dend)
							end1=r1.getTotalLength()-dstart;
					}
				}
			
			start1 = Math.max(start1, 0);
			end1 = Math.max(end1, 0);
			start2 = Math.max(start2, 0);
			end2 = Math.max(end2, 0);
			
			// correct invalid situations, e.g. when after=before=0, and start1>end2
			end2=Math.max(end2,start1-before);
			end1=Math.max(end1,start2-after);
			
			
			GenomicRegion t1 = r1.map(start1, r1.getTotalLength()-end1); // is empty if length would be negative!
			GenomicRegion t2 = r2.map(end2, r2.getTotalLength()-start2);

			GenomicRegion treg = t1.union(t2); // this is still in the induced coords system of t!
			if (treg.isEmpty()) 
				continue; // everything is clipped, might happen for very short fragments!
			AlignedReadsDataFactory fac = map.computeIfAbsent(t.map(treg), x-> 
						new AlignedReadsDataFactory(t.getData().getNumConditions(),t.getData().hasNonzeroInformation()).start());
			int tbefore,toverlap,tafter;
			if (t1.isEmpty()) {
				tbefore = toverlap = 0;
				tafter= treg.getTotalLength();
				fac.add(t.getData(), d, v->v.isFromSecondRead()?AlignedReadTrimmer.transformVariation(v,t2.getStart(),t2.getEnd(),treg):null,false);
			}
			else if (t2.isEmpty()) {
				tafter = toverlap = 0;
				tbefore = treg.getTotalLength();
				fac.add(t.getData(), d, v->v.isFromSecondRead()?null:AlignedReadTrimmer.transformVariation(v,t1.getStart(),t1.getEnd(),treg),false);
			} else {
				tbefore = treg.induce(t2.getStart());
				toverlap = treg.induce(t1.getStop())+1-tbefore;
				tafter = treg.getTotalLength()-tbefore-toverlap;
				fac.add(t.getData(), d, v->v.isFromSecondRead()?AlignedReadTrimmer.transformVariation(v,t2.getStart(),t2.getEnd(),treg):AlignedReadTrimmer.transformVariation(v,t1.getStart(),t1.getEnd(),treg),false);
			}
			
			fac.setGeometry(tbefore, toverlap, tafter);
			
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
