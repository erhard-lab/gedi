package gedi.core.data.reads;

import gedi.util.StringUtils;

public class AlignedReadsInsertion implements AlignedReadsVariation {
	
	
	private int position;
	private CharSequence sequence;
	private boolean isFromSecondRead;
	
	public AlignedReadsInsertion(int position, CharSequence sequence, boolean isFromSecondRead) {
		this.position = position;
		this.sequence = sequence;
		this.isFromSecondRead = isFromSecondRead;
	}

	@Override
	public int getPosition() {
		return position;
	}

	@Override
	public boolean isFromSecondRead() {
		return isFromSecondRead;
	}
	
	@Override
	public AlignedReadsInsertion asInsertion() {
		return this;
	}
	
	@Override
	public boolean isMismatch() {
		return false;
	}

	@Override
	public boolean isDeletion() {
		return false;
	}

	@Override
	public boolean isInsertion() {
		return true;
	}
	

	@Override
	public boolean isSoftclip() {
		return false;
	}

	@Override
	public CharSequence getReferenceSequence() {
		return "";
	}

	@Override
	public CharSequence getReadSequence() {
		return sequence;
	}
	
	@Override
	public String toString() {
		return "I"+position+sequence+(isFromSecondRead?"r":"");
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AlignedReadsInsertion other = (AlignedReadsInsertion) obj;
		if (position != other.position)
			return false;
		if (isFromSecondRead != other.isFromSecondRead)
			return false;
		if (sequence == null) {
			if (other.sequence != null)
				return false;
		} else if (!StringUtils.charsEqual(sequence,other.sequence))
			return false;
		return true;
	}

	@Override
	public int getType() {
		return DefaultAlignedReadsData.TYPE_INSERTION;
	}
	
	@Override
	public int hashCode() {
		return hashCode2();
	}

	@Override
	public AlignedReadsInsertion reposition(int newPos) {
		return new AlignedReadsInsertion(newPos, sequence, isFromSecondRead);
	}
}
