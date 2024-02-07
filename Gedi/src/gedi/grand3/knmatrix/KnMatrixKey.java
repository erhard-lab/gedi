package gedi.grand3.knmatrix;

public class KnMatrixKey {

	private String condition;
	private int subread;
	private String label;
	public KnMatrixKey(String condition, int subread, String label) {
		super();
		this.condition = condition;
		this.subread = subread;
		this.label = label;
	}
	public String getCondition() {
		return condition;
	}
	public int getSubread() {
		return subread;
	}
	public String getLabel() {
		return label;
	}
	@Override
	public String toString() {
		return "KnMatrixKey [condition=" + condition + ", subread=" + subread + ", label="
				+ label + "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((condition == null) ? 0 : condition.hashCode());
		result = prime * result + ((label == null) ? 0 : label.hashCode());
		result = prime * result + subread;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		KnMatrixKey other = (KnMatrixKey) obj;
		if (condition == null) {
			if (other.condition != null)
				return false;
		} else if (!condition.equals(other.condition))
			return false;
		if (label == null) {
			if (other.label != null)
				return false;
		} else if (!label.equals(other.label))
			return false;
		if (subread != other.subread)
			return false;
		return true;
	}
	
	
}
