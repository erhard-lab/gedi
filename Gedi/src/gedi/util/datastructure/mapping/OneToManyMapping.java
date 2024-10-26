package gedi.util.datastructure.mapping;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Set;

import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;

public class OneToManyMapping<F,T> {

	
	private HashMap<F,ArrayList<T>> map;
	private HashMap<T,F> inverse;
	
	
	public OneToManyMapping() {
	}
	
	public ArrayList<T> map(F from) {
		makeSureMap();
		return map.get(from);
	}
	public ArrayList<T> mapOrDefault(F from, ArrayList<T> defaultValue) {
		makeSureMap();
		ArrayList<T> re = map.get(from);
		if (re==null) re = defaultValue;
		return re;
	}
	public ArrayList<T> mapOrDefault(F from, T defaultValue) {
		makeSureMap();
		ArrayList<T> re = map.get(from);
		if (re==null) {
			re = new ArrayList<T>();
			re.add(defaultValue);
		}
		return re;
	}
	
	public F inverse(T to) {
		makeSureInverse();
		return inverse.get(to);
	}
	public F inverseOrDefault(T to, F defaultValue) {
		makeSureInverse();
		F re = inverse.get(to);
		if (re==null) re = defaultValue;
		return re;
	}
	
	
	public Set<F> getFromUniverse() {
		makeSureMap();
		return map.keySet();
	}
	
	public Set<T> getToUniverse() {
		makeSureInverse();
		return inverse.keySet();
	}
	
	
	private void makeSureInverse() {
		if (inverse!=null) return;
		inverse = new HashMap<T, F>();
		for (F f : map.keySet()) {
			for (T t : map.get(f))
				inverse.put(t, f);
		}
	}
	private void makeSureMap() {
		if (map!=null) return;
		map = new HashMap<F, ArrayList<T>>();
		
		for (T t : inverse.keySet()) {
			map.computeIfAbsent(inverse.get(t), x->new ArrayList<T>()).add(t);
		}
	}
	
	
	public static OneToManyMapping<String, String> fromFile(String path, String fromColumn, String toColumn) throws IOException {
		HeaderLine h = new HeaderLine();
		OneToManyMapping<String, String> re = new OneToManyMapping<String, String>();
		re.inverse = EI.lines(path).header(h).split('\t').index(a->a[h.apply(toColumn)], a->a[h.apply(fromColumn)]);
		return re;
	}
	
	
}
