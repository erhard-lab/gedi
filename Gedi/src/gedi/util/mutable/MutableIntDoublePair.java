package gedi.util.mutable;

import gedi.util.GeneralUtils;


public class MutableIntDoublePair implements Comparable<MutableIntDoublePair>, Mutable {
	
	public int Item1;
	public double Item2;
	
	public MutableIntDoublePair() {
		this(0,0);
	}
	public MutableIntDoublePair(int i1, double i2) {
		Item1 = i1;
		Item2 = i2;
	}
	
	public MutableIntDoublePair set(int i1, double i2) {
		Item1 = i1;
		Item2 = i2;
		return this;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		sb.append(String.valueOf(Item1)+",");
		sb.append(String.valueOf(Item2));
		sb.append("]");
		return sb.toString();
	}
	
	@Override
	public int hashCode() {
		int re = 0;
		re+=Integer.hashCode(Item1);
		re=31*re+Double.hashCode(Item2);
		return re;
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public boolean equals(Object obj) {
		if (obj==this) return true;
		else if (obj instanceof MutableIntDoublePair) {
			MutableIntDoublePair p = (MutableIntDoublePair) obj;
			boolean re = true;
			re &= p.Item1==Item1;
			re &= GeneralUtils.isEqual(p.Item2,Item2);
			return re;
		}
		else return false;
	}
	
	@SuppressWarnings("unchecked")
	@Override
	public int compareTo(MutableIntDoublePair o) {
		int re = 0;
		re = Integer.compare(Item1, o.Item1);
		if (re!=0) return re;
		re = Double.compare(Item2, o.Item2);
		return re;
	}
	
	public MutableIntDoublePair clone() {
		return new MutableIntDoublePair(Item1,Item2);
	}
	@Override
	public int size() {
		return 2;
	}
	@SuppressWarnings("unchecked")
	@Override
	public Object get(int index) {
		if (index==0) return Item1;
		if (index==1) return Item2;
		throw new IndexOutOfBoundsException();
	}
	@Override
	public Object set(int index, Object ele) {
		if (index==0) {
			Object re = Item1;
			Item1 = (int) ele;
			return re;
		}
		if (index==1) {
			Double re = Item2;
			Item2 = (double) ele;
			return re;
		}
		throw new IndexOutOfBoundsException();
	}
	
	public int getInt() { return Item1;}
	public double getDouble() {return Item2;}
	public int setInt(int i) { int re = Item1; Item1 = i; return re;}
	public double setDouble(double d) { double re = Item2; Item2 = d; return re;}
	
	public void incrementInt() {
		Item1++;
	}
	
}
