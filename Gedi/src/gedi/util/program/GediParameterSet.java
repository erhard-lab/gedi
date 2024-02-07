package gedi.util.program;

import java.util.ArrayList;
import java.util.HashMap;

public class GediParameterSet {
	
	private HashMap<String,GediParameter<?>> map = new HashMap<>();
	
	public void add(GediParameter<?> p) {
		map.put(p.getName(), p);
	}

	public <T> GediParameter<T> get(String name) {
		return (GediParameter<T>) map.get(name);
	}

}
