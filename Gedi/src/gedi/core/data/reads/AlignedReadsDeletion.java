package gedi.core.data.reads;

import gedi.util.StringUtils;

public class AlignedReadsDeletion implements AlignedReadsVariation{
	
	private int position;
	private CharSequence sequence;
	private boolean isFromSecondRead;
	
	public AlignedReadsDeletion(int position, CharSequence sequence, boolean isFromSecondRead) {
		this.position = position;
		this.sequence = sequence;
		this.isFromSecondRead = isFromSecondRead;
	}
	
	@Override
	public int getType() {
		return DefaultAlignedReadsData.TYPE_DELETION;
	}

	@Override
	public int getPosition() {
		return position;
	}
	
	@Override
	public AlignedReadsDeletion asDeletion() {
		return this;
	}

	@Override
	public boolean isFromSecondRead() {
		return isFromSecondRead;
	}
	
	@Override
	public boolean isMismatch() {
		return false;
	}

	@Override
	public boolean isDeletion() {
		return true;
	}

	@Override
	public boolean isInsertion() {
		return false;
	}

	@Override
	public boolean isSoftclip() {
		return false;
	}
	
	
	@Override
	public CharSequence getReferenceSequence() {
		return sequence;
	}

	@Override
	public CharSequence getReadSequence() {
		return "";
	}
	
	@Override
	public String toString() {
		return "D"+position+sequence+(isFromSecondRead?"r":"");
	}
	
	
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AlignedReadsDeletion other = (AlignedReadsDeletion) obj;
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
	public int hashCode() {
		return hashCode2();
	}

	@Override
	public AlignedReadsDeletion reposition(int newPos) {
		return new AlignedReadsDeletion(newPos, sequence, isFromSecondRead);
	}

	public boolean contains(int position) {
		return position>=this.position && position<this.position+sequence.length();
	}

}
