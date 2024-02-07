package gedi.util.oml.cps;

import gedi.util.dynamic.DynamicObject;
import gedi.util.oml.OmlInterceptor;
import gedi.util.oml.OmlNode;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.function.Function;

public class CpsList implements Function<CpsKey,DynamicObject>, OmlInterceptor {

	private LinkedHashMap<CpsKey,DynamicObject> map = new LinkedHashMap<CpsKey, DynamicObject>();
	

	void add(CpsKey key, DynamicObject obj) {
		map.put(key,obj);
	}

	
	public DynamicObject getForId(String id) {
		return apply(new CpsKey(id, null, null));
	}
	
	public DynamicObject getForClasses(HashSet<String> classes) {
		return apply(new CpsKey(null, classes, null));
	}
	
	public DynamicObject getForClasses(String... classes) {
		return apply(new CpsKey(null, new HashSet<>(Arrays.asList(classes)), null));
	}
	
	public DynamicObject getForJavaClass(Class<?> cls) {
		return apply(new CpsKey(null, null, cls));
	}
	
	public DynamicObject getForJavaClass(Object o) {
		return apply(new CpsKey(null, null, o.getClass()));
	}
	
	@Override
	public DynamicObject apply(CpsKey key) {
		ArrayList<DynamicObject> re = new ArrayList<DynamicObject>();
		for (CpsKey k : map.keySet())
			if (k.matches(key))
				re.add(map.get(k));
		return DynamicObject.cascade(re);
	}

	private int size() {
		return map.size();
	}

	@Override
	public void setObject(OmlNode node, Object o,String id, String[] classes,
			HashMap<String, Object> context) {
		if (o==null) return;
		
		DynamicObject obj = apply(new CpsKey(node.getId(), new HashSet<String>(Arrays.asList(node.getClasses())), o.getClass()));
		if (obj!=null && !obj.isNull())
			obj.applyTo(o);
		
	}
	
}
