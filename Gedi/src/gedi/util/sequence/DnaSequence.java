package gedi.util.sequence;


import cern.colt.bitvector.BitVector;
import gedi.util.SequenceUtils;


public class DnaSequence extends BitVector implements CharSequence {

	
	private static final long serialVersionUID = -2998934882733151866L;

	public DnaSequence(int len) {
		super(len*2);
	}
	
	public DnaSequence(CharSequence sequence) {
		super(countACGT(sequence)*2);
		int p = 0;
		for (int i=0; i<sequence.length(); i++) {
			char c = sequence.charAt(i);
			if (SequenceUtils.inv_nucleotides[c]<4) {
				putQuick(p*2, (SequenceUtils.inv_nucleotides[c]&1)==1);
				putQuick(p*2+1, ((SequenceUtils.inv_nucleotides[c]>>1)&1)==1);
				p++;
			}
		}
	}
	
	private static int countACGT(CharSequence sequence) {
		int re = 0;
		for (int i=0; i<sequence.length(); i++)
			if (SequenceUtils.inv_nucleotides[sequence.charAt(i)]<4)
				re++;
		return re;
	}

	public DnaSequence(char... sequence) {
		super(countACGT(sequence)*2);
		int p = 0;
		for (int i=0; i<sequence.length; i++) {
			char c = sequence[i];
			if (SequenceUtils.inv_nucleotides[c]<4) {
				putQuick(p*2, (SequenceUtils.inv_nucleotides[c]&1)==1);
				putQuick(p*2+1, ((SequenceUtils.inv_nucleotides[c]>>1)&1)==1);
				p++;
			}
		}
	}

	private static int countACGT(char[] sequence) {
		int re = 0;
		for (char c : sequence)
			if (SequenceUtils.inv_nucleotides[c]<4)
				re++;
		return re;
	}

	
	public String substring(int s, int e) {
		char[] re = new char[e-s];
		for (int i=s; i<e; i++)
			re[i-s] = charAt(i);
		return new String(re);
	}
	
	@Override
	public String toString() {
		return getSequence();
	}
	
	public int getSequenceLength() {
		return size()/2;
	}
	
	
	public String getSequence() {
		char[] re = new char[getSequenceLength()];
		for (int i=0; i<re.length; i++)
			re[i] = charAt(i);
		return new String(re);
	}
	

	@Override
	public char charAt(int index) {
		int r=0;
		if (getQuick(index*2))
			r|=1;
		if (getQuick(index*2+1))
			r|=2;
		return SequenceUtils.nucleotides[r];
	}

	@Override
	public int length() {
		return getSequenceLength();
	}

	@Override
	public DnaSequence subSequence(int start, int end) {
		return new DnaSequence(getSequence().substring(start, end));
	}
	
	@Override
	public boolean equals(Object obj) {
		if (obj==this) return true;
		if (obj==null) return false;
		if (!(obj instanceof CharSequence)) return false;
		CharSequence s = (CharSequence) obj;
		if (s.length()!=length()) return false;
		for (int i=0; i<s.length(); i++) 
			if (s.charAt(i)!=charAt(i)) return false;
		return true;
	}
	
}
