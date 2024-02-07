package gedi.util.datastructure.charsequence;

import java.util.Arrays;

import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import cern.colt.bitvector.BitVector;

public class CharArrayCharSequence implements CharSequence {

	private char[] seq;
	
	public CharArrayCharSequence(char[] seq) {
		this.seq = seq;
	}

	@Override
	public int length() {
		return seq.length;
	}

	@Override
	public char charAt(int index) {
		return seq[index];
	}
	

	@Override
	public CharSequence subSequence(int start, int end) {
		return new CharArrayCharSequence(Arrays.copyOfRange(seq, start, end));
	}
	
	@Override
	public String toString() {
		return StringUtils.toString(this);
	}
	
	@Override
	public boolean equals(Object obj) {
		return StringUtils.equals(this, obj);
	}
	
	@Override
	public int hashCode() {
		return StringUtils.hashCode(this);
	}

	
}
