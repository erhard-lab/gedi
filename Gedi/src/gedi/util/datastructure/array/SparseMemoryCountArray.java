package gedi.util.datastructure.array;

import gedi.util.functions.BiIntConsumer;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;

import java.io.IOException;

/**
 * serialized as CInt!
 * @author flo
 *
 */
public class SparseMemoryCountArray extends MemoryIntegerArray {
	
	
	public SparseMemoryCountArray() {
	}
	
	public SparseMemoryCountArray(int[] wrap) {
		super(wrap);
	}
	public SparseMemoryCountArray(int length) {
		super(length);
	}

	
	public void forEachNonZero(BiIntConsumer action) {
		for (int i=0; i<length(); i++)
			if (!isZero(i))
				action.accept(i, getInt(i));
	}
	
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		int zeros = 0;
		for (int i=0; i<length(); i++)
			if (getInt(i)==0)
				zeros++;
		if (2*length()>(length()-zeros)*5+1) {
			// sparse
			out.putInt(-length());
			out.putCInt(zeros);
			for (int i=0; i<length(); i++)
				if (getInt(i)!=0) {
					out.putCInt(i);
					out.putCInt(getInt(i));
				}
		} else {
			out.putInt(length());
			for (int i=0; i<length(); i++)
				out.putCInt(getInt(i));
		}
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		int len = in.getInt();
		if (len>=0) {
			if (a==null || a.length!=len)
				a = new int[len];
			for (int i=0; i<a.length; i++)
				a[i] = in.getCInt();
		} else {
			len = -len;
			if (a==null || a.length!=len)
				a = new int[len];
			int elems = a.length-in.getCInt();
			for (int i=0; i<elems; i++) {
				int ind = in.getCInt();
				a[ind] = in.getCInt();
			}
		}
	}
	
}
