package gedi.util.datastructure.tree.redblacktree;

public class UnitInterval implements Interval{

	private int pos;
	
	public UnitInterval(int pos) {
		this.pos = pos;
	}

	@Override
	public int getStart() {
		return pos;
	}

	@Override
	public int getStop() {
		return pos;
	}

	@Override
	public int hashCode() {
		return pos;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		UnitInterval other = (UnitInterval) obj;
		if (pos != other.pos)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return pos+"";
	}

	
	

}
