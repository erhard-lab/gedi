package gedi.gui.gtracks.style;

import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;

import gedi.app.classpath.ClassPathCache;
import gedi.util.StringUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;

public class GTrackStyles {
	
	private LinkedHashMap<String, Object>[] styles;
	private HashMap<String,LinkedHashMap<String, Object>> nameMap;
	
	public GTrackStyles(DynamicObject dob) {
		styles = EI.wrap(dob.asArray()).map(d->d.asMap(GTrackStyles::read)).toArray((Class)LinkedHashMap.class);
		nameMap = EI.wrap(styles).index(m->m.get("name").toString());
	}
	
	
	public <T> T getStyle(int index, String style) {
		return (T) styles[index].get(style);
	}
	
	public <T> T getStyle(String name, String style) {
		return (T) nameMap.get(name).get(style);
	}

	
	public <T> T getStyle(int index, String style, T defaultVal) {
		return (T) styles[index].getOrDefault(style,defaultVal);
	}
	
	public <T> T getStyle(String name, String style, T defaultVal) {
		return (T) nameMap.get(name).getOrDefault(style,defaultVal);
	}

	
	private static Map<String,GTrackStyleParser> parsers = Collections.synchronizedMap(new HashMap<String,GTrackStyleParser>());
	
	private static GTrackStyleParser obtainParser(String name) {
		try {
			return (GTrackStyleParser) ClassPathCache.getInstance().getClass("GTrack"+StringUtils.getCamelCase(name,true)+"StyleParser").newInstance();
		} catch (InstantiationException | IllegalAccessException e) {
			throw new RuntimeException("Cannot set GTrack style! Implement GTrack"+StringUtils.getCamelCase(name,true)+"StyleParser!",e);
		}
	}
	
	public static Object read(String name, DynamicObject value) {
		GTrackStyleParser parser = parsers.computeIfAbsent(name, GTrackStyles::obtainParser);
		return parser.apply(value);
	}

}

