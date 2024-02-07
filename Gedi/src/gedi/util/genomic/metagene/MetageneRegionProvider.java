package gedi.util.genomic.metagene;

import java.util.function.Function;

import gedi.core.data.annotation.Transcript;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.ImmutableReferenceGenomicRegion;

public interface MetageneRegionProvider {
	
	

	public static GenomicRegion getRegion(ImmutableReferenceGenomicRegion<?> g, GenomicRegionPosition pos1, int offset1, GenomicRegionPosition pos2, int offset2) {
		int p1 = pos1.position(g,offset1);
		int p2 = pos2.position(g,offset2);
		if (g.getReference().isMinus()) return new ArrayGenomicRegion(p2,p1);
		return new ArrayGenomicRegion(p1,p2);
	}
	
	public static MetageneRegionProvider make(Function<ImmutableReferenceGenomicRegion<?>,GenomicRegion> reg, String unit, String start, String stop) {
		return new MetageneRegionProvider() {
			@Override
			public String getUnit() {
				return unit;
			}
			@Override
			public String getStop() {
				return stop;
			}
			@Override
			public String getStart() {
				return start;
			}
			@Override
			public GenomicRegion getRegion(ImmutableReferenceGenomicRegion<?> g) {
				return reg.apply(g);
			}
		};
	}
	
	public static MetageneRegionProvider range(GenomicRegionPosition pos1, int offset1, GenomicRegionPosition pos2, int offset2, String unit, String start, String stop) {
		return new MetageneRegionProvider() {
			@Override
			public String getUnit() {
				return unit;
			}
			@Override
			public String getStop() {
				return stop;
			}
			@Override
			public String getStart() {
				return start;
			}
			@Override
			public GenomicRegion getRegion(ImmutableReferenceGenomicRegion<?> g) {
				return MetageneRegionProvider.getRegion(g,pos1,offset1,pos2,offset2);
			}
		};
	}
	
	public static MetageneRegionProvider full() {
		return range(GenomicRegionPosition.FivePrime,0,GenomicRegionPosition.ThreePrime,0,"%","0","100");
	}

	public static MetageneRegionProvider inner(int removePrefix, int removeSuffix) {
		return range(GenomicRegionPosition.FivePrime,removePrefix,GenomicRegionPosition.ThreePrime,-removeSuffix,"%",removePrefix+"","-"+removeSuffix);
	}

	public static MetageneRegionProvider upstream(int len) {
		return range(GenomicRegionPosition.FivePrime,-1-len,GenomicRegionPosition.FivePrime,-1,"bp","-"+len,"0");
	}
	
	public static MetageneRegionProvider prefix(int len) {
		return range(GenomicRegionPosition.FivePrime,0,GenomicRegionPosition.FivePrime,len,"bp","0",""+len);
	}
	
	public static MetageneRegionProvider suffix(int len) {
		return range(GenomicRegionPosition.ThreePrime,-len,GenomicRegionPosition.ThreePrime,0,"bp","-"+len,"0");
	}
	
	public static MetageneRegionProvider downstream(int len) {
		return range(GenomicRegionPosition.ThreePrime,1,GenomicRegionPosition.ThreePrime,1+len,"bp","0",""+len);
	}
	
	public static MetageneRegionProvider utr5() {
		return make(g->((Transcript) g.getData()).get5Utr(g).getRegion(),"%","0","100");
	}
	
	public static MetageneRegionProvider cds() {
		return make(g->((Transcript) g.getData()).getCds(g).getRegion(),"%","0","100");
	}
	
	public static MetageneRegionProvider utr3() {
		return make(g->((Transcript) g.getData()).get3Utr(g).getRegion(),"%","0","100");
	}
	
	
	public static MetageneRegionProvider innerCds(int removePrefix, int removeSuffix) {
		return range(GenomicRegionPosition.Start,removePrefix,GenomicRegionPosition.Stop,-removeSuffix,"%",removePrefix+"","-"+removeSuffix);
	}

	public static MetageneRegionProvider upstreamCds(int len) {
		return range(GenomicRegionPosition.Start,-1-len,GenomicRegionPosition.Start,-1,"bp","-"+len,"0");
	}
	
	public static MetageneRegionProvider prefixCds(int len) {
		return range(GenomicRegionPosition.Start,0,GenomicRegionPosition.Start,len,"bp","0",""+len);
	}
	
	public static MetageneRegionProvider suffixCds(int len) {
		return range(GenomicRegionPosition.Stop,-len,GenomicRegionPosition.Stop,0,"bp","-"+len,"0");
	}
	
	public static MetageneRegionProvider downstreamCds(int len) {
		return range(GenomicRegionPosition.Stop,1,GenomicRegionPosition.Stop,1+len,"bp","0",""+len);
	}
	
	default ReferenceSequence getReference(ImmutableReferenceGenomicRegion<?> g) {
		return g.getReference();
	}
	GenomicRegion getRegion(ImmutableReferenceGenomicRegion<?> g);
	
	String getUnit();
	String getStart();
	String getStop();
	
}
