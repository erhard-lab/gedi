package gedi.grand3.knmatrix;

import java.util.Arrays;

public final class SubreadKNKey {
	
	public int[] k;
	public int[] n;
	
	public SubreadKNKey(int numSubreads) {
		this.k = new int[numSubreads];
		this.n = new int[numSubreads];
	}
	
	private SubreadKNKey() {}

	public SubreadKNKey set(int[] k, int[] n) {
		this.k = k;
		this.n = n;
		return this;
	}
	
	public SubreadKNKey clone() {
		SubreadKNKey re = new SubreadKNKey();
		re.k = k.clone();
		re.n = n.clone();
		return re;
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(k);
		result = prime * result + Arrays.hashCode(n);
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
		SubreadKNKey other = (SubreadKNKey) obj;
		if (!Arrays.equals(k, other.k))
			return false;
		if (!Arrays.equals(n, other.n))
			return false;
		return true;
	}
	
	
	
	
}