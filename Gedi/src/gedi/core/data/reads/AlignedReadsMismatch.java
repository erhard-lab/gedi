package gedi.core.data.reads;

public class AlignedReadsMismatch implements AlignedReadsVariation {

	private int position;
	private CharSequence genomic;
	private CharSequence read;
	private boolean isFromSecondRead;
	
	public AlignedReadsMismatch(int position, CharSequence genomic, CharSequence read, boolean isFromSecondRead) {
		this.position = position;
		this.genomic = genomic;
		this.read = read;
		this.isFromSecondRead = isFromSecondRead;
	}

	@Override
	public int getPosition() {
		return position;
	}

	@Override
	public int getType() {
		return DefaultAlignedReadsData.TYPE_MISMATCH;
	}
	
	@Override
	public AlignedReadsMismatch asMismatch() {
		return this;
	}
	
	@Override
	public boolean isMismatch() {
		return true;
	}

	@Override
	public boolean isDeletion() {
		return false;
	}

	@Override
	public boolean isInsertion() {
		return false;
	}
	
	@Override
	public boolean isFromSecondRead() {
		return isFromSecondRead;
	}


	@Override
	public boolean isSoftclip() {
		return false;
	}
	
	
	@Override
	public CharSequence getReferenceSequence() {
		return genomic;
	}

	@Override
	public CharSequence getReadSequence() {
		return read;
	}
	
	@Override
	public String toString() {
		return "M"+position+genomic+read+(isFromSecondRead?"r":"");
	}

	@Override
	public int hashCode() {
		return hashCode2();
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		AlignedReadsMismatch other = (AlignedReadsMismatch) obj;
		if (genomic == null) {
			if (other.genomic != null)
				return false;
		} else if (!genomic.equals(other.genomic))
			return false;
		if (isFromSecondRead != other.isFromSecondRead)
			return false;
		if (position != other.position)
			return false;
		if (read == null) {
			if (other.read != null)
				return false;
		} else if (!read.equals(other.read))
			return false;
		return true;
	}

	@Override
	public AlignedReadsMismatch reposition(int newPos) {
		return new AlignedReadsMismatch(newPos, genomic, read, isFromSecondRead);
	}

	
	
}
