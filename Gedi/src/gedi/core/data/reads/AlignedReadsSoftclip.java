package gedi.core.data.reads;

import gedi.util.StringUtils;

public class AlignedReadsSoftclip implements AlignedReadsVariation {

	private boolean p5;
	private CharSequence read;
	private boolean isFromSecondRead;
	
	public AlignedReadsSoftclip(boolean p5, CharSequence read, boolean isFromSecondRead) {
		this.p5 = p5;
		this.read = read;
		this.isFromSecondRead = isFromSecondRead;
	}

	@Override
	public int getPosition() {
		return p5?0:Integer.MAX_VALUE;
	}

	@Override
	public int getType() {
		return DefaultAlignedReadsData.TYPE_SOFTCLIP;
	}
	
	@Override
	public AlignedReadsSoftclip asSoftclip() {
		return this;
	}
	
	@Override
	public boolean isMismatch() {
		return false;
	}

	@Override
	public boolean isFromSecondRead() {
		return isFromSecondRead;
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
	public boolean isSoftclip() {
		return true;
	}
	
	
	@Override
	public CharSequence getReferenceSequence() {
		return "";
	}

	@Override
	public CharSequence getReadSequence() {
		return read;
	}
	
	@Override
	public String toString() {
		return (p5?"5p":"3p")+read+(isFromSecondRead?"r":"");
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
		AlignedReadsSoftclip other = (AlignedReadsSoftclip) obj;
		if (p5 != other.p5)
			return false;
		if (isFromSecondRead != other.isFromSecondRead)
			return false;
		if (read == null) {
			if (other.read != null)
				return false;
		} else if (!StringUtils.charsEqual(read,other.read))
			return false;
		return true;
	}

	@Override
	public AlignedReadsVariation reposition(int newPos) {
		return this;
	}
	
	
}
