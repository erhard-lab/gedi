package gedi.util.datastructure.array;

import gedi.util.io.randomaccess.BinaryReader;

import java.io.IOException;
import java.util.Arrays;

public class MemoryDoubleArray extends DoubleArray {
	
	private double[] a;
	
	public MemoryDoubleArray() {
	}
	
	public MemoryDoubleArray(double[] wrap) {
		this.a = wrap;
	}
	public MemoryDoubleArray(int length) {
		if (length>=0) a = new double[length];
	}

	@Override
	public void setDouble(int index, double value) {
		a[index] = value;
	}

	@Override
	public boolean isReadOnly() {
		return false;
	}
	
	@Override
	public NumericArray clear() {
		Arrays.fill(a, 0);
		return this;
	}

	@Override
	public int length() {
		return a.length;
	}

	@Override
	public double getDouble(int index) {
		return a[index];
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		int len = in.getInt();
		if (a==null || a.length!=len)
			a = new double[len];
		for (int i=0; i<a.length; i++)
			a[i] = in.getDouble();
	}
	
	
	@Override
	public int compare(int index1, NumericArray a2, int index2) {
		return Double.compare(a[index1], a2.getDouble(index2));
	}

	@Override
	public int compare(int index1, int index2) {
		return Double.compare(a[index1], a[index2]);
	}
	
	@Override
	public int compareInCum(int index1, int index2) {
		double b1 = index1==0?a[index1]:(a[index1]-a[index1-1]);
		double b2 = index2==0?a[index2]:(a[index2]-a[index2-1]);
		return Double.compare(b1,b2);
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
        	double tmp = a[from];
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
		double tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
		return this;
	}

	@Override
	public NumericArray copy(NumericArray fromArray, int fromIndex, int to) {
		a[to]=fromArray.getDouble(fromIndex);
		return this;
	}

	@Override
	public NumericArray copyRange(int start, NumericArray dest, int destOffset,
			int len) {
		if (dest instanceof MemoryDoubleArray)
			System.arraycopy(a, start, ((MemoryDoubleArray)dest).a, destOffset, len);
		else
			return super.copyRange(start, dest, destOffset, len);
		return this;
	}

}
