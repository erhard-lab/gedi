package gedi.util.datastructure.array;

import gedi.util.io.randomaccess.BinaryReader;

import java.io.IOException;
import java.util.Arrays;

public class MemoryLongArray extends LongArray {
	
	private long[] a;
	
	public MemoryLongArray() {
	}
	
	public MemoryLongArray(long[] wrap) {
		this.a = wrap;
	}
	public MemoryLongArray(int length) {
		if (length>=0) a = new long[length];
	}

	@Override
	public void setLong(int index, long value) {
		a[index] = value;
	}

	@Override
	public boolean isReadOnly() {
		return false;
	}

	@Override
	public int length() {
		return a.length;
	}
	@Override
	public NumericArray clear() {
		Arrays.fill(a, 0);
		return this;
	}

	
	@Override
	public long getLong(int index) {
		return a[index];
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		int len = in.getInt();
		if (a==null || a.length!=len)
			a = new long[len];
		for (int i=0; i<a.length; i++)
			a[i] = in.getLong();
	}
	
	
	@Override
	public int compare(int index1, NumericArray a2, int index2) {
		return Long.compare(a[index1], a2.getLong(index2));
	}

	@Override
	public int compare(int index1, int index2) {
		return Long.compare(a[index1], a[index2]);
	}
	
	@Override
	public int compareInCum(int index1, int index2) {
		long b1 = index1==0?a[index1]:(a[index1]-a[index1-1]);
		long b2 = index2==0?a[index2]:(a[index2]-a[index2-1]);
		return Long.compare(b1,b2);
	}
	
	
	@Override
	public NumericArray copy(int from, int to) {
		a[to] = a[from];
		return this;
	}
	
	@Override
	public NumericArray sort(int from, int to) {
		Arrays.sort(a, from, to);
		return this;
	}
	
	@Override
	public NumericArray reverse(int from, int to) {
		to--;
        while (from< to) {
        	long tmp = a[from];
    		a[from++] = a[to];
    		a[to--] = tmp;
        }
        return this;
	}
	

	@Override
	public NumericArray cumSum(int from, int to) {
		for (int i=from+1; i<to; i++)
			a[i] += a[i-1];
		return this;
	}

	@Override
	public NumericArray deCumSum(int from, int to) {
		for (int i=to-1; i>from; i--)
			a[i] -= a[i-1];
		return this;
	}

	@Override
	public NumericArray switchElements(int i, int j) {
		long tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
		return this;
	}

	@Override
	public NumericArray copy(NumericArray fromArray, int fromIndex, int to) {
		a[to]=fromArray.getLong(fromIndex);
		return this;
	}

	@Override
	public NumericArray copyRange(int start, NumericArray dest, int destOffset,
			int len) {
		if (dest instanceof MemoryLongArray)
			System.arraycopy(a, start, ((MemoryLongArray)dest).a, destOffset, len);
		else
			return super.copyRange(start, dest, destOffset, len);
		return this;
	}

}
