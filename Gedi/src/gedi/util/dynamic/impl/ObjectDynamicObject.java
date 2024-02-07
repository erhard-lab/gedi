package gedi.util.dynamic.impl;

import gedi.util.ReflectionUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.orm.Orm;
import gedi.util.orm.Orm.OrmInfo;

import java.util.Arrays;
import java.util.Collection;



public class ObjectDynamicObject implements DynamicObject {
	

	private Object o;
	private OrmInfo orm;
	
	public ObjectDynamicObject(Object o) {
		this.o = o;
		this.orm = Orm.getInfo(o.getClass());
	}

	@Override
	public Class<?> getType() {
		return Object.class;
	}


	@Override
	public boolean hasProperty(String property) {
		return orm.getIndex(property)!=-1;
	}
	
	@Override
	public String asString() {
		return "";
	}

	@Override
	public int asInt() {
		return 0;
	}

	@Override
	public double asDouble() {
		return 0;
	}

	@Override
	public boolean asBoolean() {
		return false;
	}


	@Override
	public String toString() {
		return o.toString();
	}


	@Override
	public DynamicObject getEntry(int index) {
		return DynamicObject.getEmpty();
	}

	@Override
	public DynamicObject getEntry(String property) {
		try {
			return DynamicObject.from((Object)ReflectionUtils.get(o, property));
		} catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
			throw new RuntimeException("Cannot wrap object in DynamicObject!",e);
		}
	}

	@Override
	public Collection<String> getProperties() {
		return Arrays.asList(orm.getDeclaredNames());
	}

	@Override
	public int length() {
		return 0;
	}
	
}
