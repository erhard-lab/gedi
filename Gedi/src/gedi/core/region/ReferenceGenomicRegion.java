package gedi.core.region;

import java.util.function.UnaryOperator;

import gedi.core.data.annotation.GenomicRegionMappable;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.util.StringUtils;



/**
 * Contracts:
 * compareTo does not consider the data
 * 
 * equals and hashcode do!
 * 
 * @author erhard
 *
 * @param <D>
 */
public interface ReferenceGenomicRegion<D> extends Comparable<ReferenceGenomicRegion<D>>  {

	ReferenceSequence getReference();
	GenomicRegion getRegion();
	D getData();
	
	default MutableReferenceGenomicRegion<D> toMutable() {
		return new MutableReferenceGenomicRegion<D>().set(getReference(), getRegion(), getData());
	}
	
	default ImmutableReferenceGenomicRegion<D> toImmutable() {
		return new ImmutableReferenceGenomicRegion<D>(getReference(), getRegion(), getData());
	}

	default ImmutableReferenceGenomicRegion<D> toImmutable(UnaryOperator<D> dataToImmutable) {
		return new ImmutableReferenceGenomicRegion<D>(getReference(), getRegion(), dataToImmutable.apply(getData()));
	}
	
	default <T> ImmutableReferenceGenomicRegion<T> toImmutable(T data) {
		return new ImmutableReferenceGenomicRegion<T>(getReference(), getRegion(), data);
	}
	
	boolean isMutable();
	
	default int compareTo(ReferenceGenomicRegion<D> o) {
		int re = getReference().compareTo(o.getReference());
		if (re==0)
			re = getRegion().compareTo(o.getRegion());
		return re;
	}
	
	default int compareToIgnoreStrand(ReferenceGenomicRegion<D> o) {
		int re = ReferenceSequence.compareChromosomeNames(getReference().getName(),o.getReference().getName());
		if (re==0)
			re = getRegion().compareTo(o.getRegion());
		return re;
	}
	
	default boolean intersects(ReferenceGenomicRegion<?> o) {
		return getReference().equals(o.getReference()) && getRegion().intersects(o.getRegion());
	}
	
	default boolean contains(ReferenceGenomicRegion<?> o) {
		return getReference().equals(o.getReference()) && getRegion().contains(o.getRegion());
	}
	
	default int map(int position) {
		return getReference().getStrand()==Strand.Minus?getRegion().map(getRegion().getTotalLength()-1-position):getRegion().map(position);
	}
	default GenomicRegion map(GenomicRegion region) {
		return getReference().getStrand()==Strand.Minus?getRegion().map(region.reverse(getRegion().getTotalLength())):getRegion().map(region);
	}
	default GenomicRegion mapMaybeOutSide(GenomicRegion region) {
		return getReference().getStrand()==Strand.Minus?getRegion().mapMaybeOutside(region.reverse(getRegion().getTotalLength())):getRegion().mapMaybeOutside(region);
	}
	default int mapMaybeOutSide(int pos) {
		return getReference().getStrand()==Strand.Minus?getRegion().mapMaybeOutside(getRegion().getTotalLength()-(pos+1)):getRegion().mapMaybeOutside(pos);
	}

	default <T> ImmutableReferenceGenomicRegion<T> map(ReferenceGenomicRegion<T> region) {
		ReferenceSequence ref = getReference();
		if (region.getReference().getStrand()==Strand.Minus) ref = ref.toOppositeStrand();
		else if (region.getReference().getStrand()==Strand.Independent) ref = ref.toStrandIndependent();
		T data = region.getData();
		if (data instanceof GenomicRegionMappable)
			data = (T) ((GenomicRegionMappable) data).map(this);
		return new ImmutableReferenceGenomicRegion<T>(ref, map(region.getRegion()),data);
	}
	
	default <T> ImmutableReferenceGenomicRegion<T> induce(ReferenceGenomicRegion<T> region, String referenceName) {
		if (!getReference().getName().equals(region.getReference().getName()))
			throw new IllegalArgumentException("Cannot induce a region from another reference!");
		
		ReferenceSequence ref = Chromosome.obtain(referenceName, Strand.Plus);
		if (region.getReference().getStrand()!=getReference().getStrand()) ref = ref.toOppositeStrand();
		else if (region.getReference().getStrand()==Strand.Independent) ref = ref.toStrandIndependent();
		
		T data = region.getData();
		if (data instanceof GenomicRegionMappable)
			data = (T) ((GenomicRegionMappable) data).induce(this);
		
		return new ImmutableReferenceGenomicRegion<T>(ref, induce(region.getRegion()),data);
	}
	
	
	default int induce(int position) {
		return getReference().getStrand()==Strand.Minus?getRegion().getTotalLength()-1-getRegion().induce(position):getRegion().induce(position);
	}
	default int induceMaybeOutside(int position) {
		return getReference().getStrand()==Strand.Minus?getRegion().getTotalLength()-1-getRegion().induceMaybeOutside(position):getRegion().induceMaybeOutside(position);
	}
	default ArrayGenomicRegion induce(GenomicRegion region) {
		return getReference().getStrand()==Strand.Minus?getRegion().induce(region).reverse(getRegion().getTotalLength()):getRegion().induce(region);
	}
	
	
	/**
	 * If the 5' end of region is downstream of the 3' end of this.
	 * @param region
	 * @return
	 */
	default boolean isDownstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getStart()>region.getStop();
		}
		return getRegion().getStop()<region.getStart();
	}

	/**
	 * If the 3' end of region is upstream of the 5' end of this.
	 * @param region
	 * @return
	 */
	default boolean isUpstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getStop()<region.getStart();
		}
		return getRegion().getStart()>region.getStop();
	}

	/**
	 * If the 5' end of region is downstream of the 5' end of this.
	 * @param region
	 * @return
	 */
	default boolean startsDownstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getEnd()>region.getEnd();
		}
		return getRegion().getStart()<region.getStart();
	}

	/**
	 * If the 5' end of region is upstream of the 5' end of this.
	 * @param region
	 * @return
	 */
	default boolean startsUpstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getEnd()<region.getEnd();
		}
		return getRegion().getStart()>region.getStart();
	}
	
	
	/**
	 * If the 3' end of region is downstream of the 3' end of this.
	 * @param region
	 * @return
	 */
	default boolean endsDownstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getStart()>region.getStart();
		}
		return getRegion().getEnd()<region.getEnd();
	}
	
	/**
	 * If the 5' end of region is upstream of the 5' end of this.
	 * @param region
	 * @return
	 */
	default boolean endsUpstream(GenomicRegion region) {
		if (getReference().getStrand()==Strand.Minus) {
			return getRegion().getStart()<region.getStart();
		}
		return getRegion().getEnd()>region.getEnd();
	}
	
	default ImmutableReferenceGenomicRegion<D> getUpstream(int len) {
		GenomicRegion reg;
		if (getReference().getStrand()==Strand.Minus) 
			reg=new ArrayGenomicRegion(getRegion().getEnd(),getRegion().getEnd()+len);
		else
			reg=new ArrayGenomicRegion(getRegion().getStart()-len,getRegion().getStart());
		return new ImmutableReferenceGenomicRegion<>(getReference(), reg,getData());
	}
	
	default ImmutableReferenceGenomicRegion<D> getDownstream(int len) {
		GenomicRegion reg;
		if (getReference().getStrand()==Strand.Minus)
			reg=new ArrayGenomicRegion(getRegion().getStart()-len,getRegion().getStart());
		else
			reg=new ArrayGenomicRegion(getRegion().getEnd(),getRegion().getEnd()+len);
			
		return new ImmutableReferenceGenomicRegion<>(getReference(), reg,getData());
	}
	
	
	default int hashCode2() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((getData() == null) ? 0 : getData().hashCode());
		result = prime * result
				+ ((getReference() == null) ? 0 : getReference().hashCode());
		result = prime * result + ((getRegion() == null) ? 0 : getRegion().hashCode());
		return result;
	}

	default boolean equals2(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof ReferenceGenomicRegion))
			return false;
		ReferenceGenomicRegion<D> other = (ReferenceGenomicRegion<D>) obj;
		if (getData() == null) {
			if (other.getData() != null)
				return false;
		} else if (!getData().equals(other.getData()))
			return false;
		if (getReference() == null) {
			if (other.getReference() != null)
				return false;
		} else if (!getReference().equals(other.getReference()))
			return false;
		if (getRegion() == null) {
			if (other.getRegion() != null)
				return false;
		} else if (!getRegion().equals(other.getRegion()))
			return false;
		return true;
	}

	default String toString2() {
		String reg = getRegion().toString();
		if (getData() instanceof AlignedReadsData) 
			reg = getRegion().toString((ReferenceGenomicRegion)this);
		return getData()==null?getReference()+":"+getRegion():getReference()+":"+reg+"\t"+StringUtils.toString(getData());
	}
	
	default String toLocationString() {
		String reg = getRegion().toRegionString();
		if (getData() instanceof AlignedReadsData) 
			reg = getRegion().toString((ReferenceGenomicRegion)this);
		return getReference().toPlusMinusString()+":"+reg;
	}
	
	default String toLocationStringRemovedIntrons() {
		return getReference().toPlusMinusString()+":"+getRegion().removeIntrons().toRegionString();
	}
	
	

}
