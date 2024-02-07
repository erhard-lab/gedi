package gedi.core.region;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.util.StringUtils;

public class ImmutableReferenceGenomicRegion<D> implements ReferenceGenomicRegion<D> {

	private ReferenceSequence reference;
	private GenomicRegion region;
	private D data;
	
	private transient int hashcode = -1;

	public ImmutableReferenceGenomicRegion(ReferenceSequence reference,
			GenomicRegion region) {
		this(reference,region,null);
	}
	
	public ImmutableReferenceGenomicRegion(ReferenceSequence reference,
			GenomicRegion region, D data) {
		this.reference = reference;
		this.region = region;
		this.data = data;
	}
	
	@Override
	public boolean isMutable() {
		return false;
	}

	public ReferenceSequence getReference() {
		return reference;
	}

	public GenomicRegion getRegion() {
		return region;
	}
	
	public D getData() {
		return data;
	}
	
	@Override
	public ImmutableReferenceGenomicRegion<D> toImmutable() {
		return this;
	}

	public ImmutableReferenceGenomicRegion<D> toStrandIndependent() {
		return new ImmutableReferenceGenomicRegion<>(getReference().toStrandIndependent(), getRegion(),getData());
	}
	
	public ImmutableReferenceGenomicRegion<D> toPlusStrand() {
		return new ImmutableReferenceGenomicRegion<>(getReference().toPlusStrand(), getRegion(),getData());
	}
	
	public ImmutableReferenceGenomicRegion<D> toOppositeStrand() {
		return new ImmutableReferenceGenomicRegion<>(getReference().toOppositeStrand(), getRegion(),getData());
	}
	
	public ImmutableReferenceGenomicRegion<D> toMinusStrand() {
		return new ImmutableReferenceGenomicRegion<>(getReference().toMinusStrand(), getRegion(),getData());
	}
	

	public ImmutableReferenceGenomicRegion<D> toStrand(Strand strand) {
		return new ImmutableReferenceGenomicRegion<>(getReference().toStrand(strand), getRegion(),getData());
	}
	
	@Override
	public String toString() {
		return toString2();
	}
	
	@Override
	public int hashCode() {
		if (hashcode==-1)
			hashcode = hashCode2();
		return hashcode;
	}

	@Override
	public boolean equals(Object obj) {
		return equals2(obj);
	}
	
	public static boolean canParse(String pos) {
		int ind = pos.indexOf('\t');
		if (ind>=0) pos.substring(0, ind);
		String[] p = StringUtils.split(pos,':');
		if (p.length!=2) return false;
		return GenomicRegion.parse(p[1])!=null;
	}

	public static <D> ImmutableReferenceGenomicRegion<D> parse(String pos) {
		int ind = pos.indexOf('\t');
		if (ind>=0)
			return parse(pos.substring(0, ind),(D) new NameAnnotation(pos.substring(ind+1)));
		return parse(pos,null);
	}
		
	public static <D> ImmutableReferenceGenomicRegion<D> parse(String pos, D data) {
		String[] p = StringUtils.split(pos,':');
		if (p.length!=2) throw new RuntimeException("Cannot parse as location: "+pos);
		return new ImmutableReferenceGenomicRegion<D>(Chromosome.obtain(p[0]), GenomicRegion.parse(p[1]),data);
	}

	
	public static <D> ImmutableReferenceGenomicRegion<D> parse(Genomic g, String pos) {
		return parse(g,pos,null);
	}
		
	public static <D> ImmutableReferenceGenomicRegion<D> parse(Genomic g, String pos, D data) {
		if (!pos.contains(":")) {
			
			if ((pos.endsWith("-") || pos.endsWith("+")) && g.getSequenceNames().contains(pos.substring(0,pos.length()-1))) {
				return new ImmutableReferenceGenomicRegion<D>(Chromosome.obtain(pos),g.getRegionOfReference(pos.substring(0,pos.length()-1)),data);
			}
			
			
			boolean opp = pos.endsWith("-");
			if (opp) pos = pos.substring(0,pos.length()-1);
			
			ImmutableReferenceGenomicRegion<?> rgr = g.getNameIndex().get(pos).toImmutable();
			if (rgr!=null) {
				rgr = rgr.toMutable().transformRegion(r->r.extendFront(1)).toImmutable();// workaround for IndexGenome bug
				if (opp) rgr=rgr.toOppositeStrand();
				pos = rgr.toLocationString();
			} else {
				rgr = g.getGeneMapping().apply(pos).toImmutable();
				if (rgr!=null) {
					if (opp) rgr=rgr.toOppositeStrand();
					pos = rgr.toLocationString();
				} else {
					rgr = g.getTranscriptMapping().apply(pos).toImmutable();
					if (rgr!=null) {
						if (opp) rgr=rgr.toOppositeStrand();
						pos = rgr.toLocationString();
					} else {
						ReferenceSequence ref = pos==null?g.getTranscripts().getReferenceSequences().iterator().next():Chromosome.obtain(pos);
						GenomicRegion reg = g.getTranscripts().getTree(ref).getRoot().getKey().removeIntrons();
						reg = reg.extendAll(reg.getTotalLength()/3, reg.getTotalLength()/3);
						pos = ref.toPlusMinusString()+":"+reg.toRegionString();
					}
				}
			}
		}
		return parse(pos,data);
	}

	
	
	
}
