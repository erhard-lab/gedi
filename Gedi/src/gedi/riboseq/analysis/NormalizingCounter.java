package gedi.riboseq.analysis;

import java.util.HashMap;
import java.util.Set;


class NormalizingCounter<T> {

	private HashMap<T,CountObject> values = new HashMap<>();
	
	public void count(T key, double val) {
		CountObject co = values.computeIfAbsent(key, x->new CountObject());
		co.counted++;
		co.value+=val;
	}
	
	public void add(NormalizingCounter<T> c) {
		for (T key : c.values.keySet()) {
			CountObject oo = c.values.get(key);
			CountObject co = values.computeIfAbsent(key, x->new CountObject());
			co.counted+=oo.counted;
			co.value+=oo.value;
		}
	}
	
	public Set<T> keySet() {
		return values.keySet();
	}
	
	public double value(T key) {
		CountObject oo = values.get(key);
		if (oo==null) return 0;
		return oo.value/oo.counted;
	}
	
	public int count(T key) {
		CountObject oo = values.get(key);
		if (oo==null) return 0;
		return oo.counted;
	}
	
	private static class CountObject {
		int counted;
		double value;
	}
	
}
