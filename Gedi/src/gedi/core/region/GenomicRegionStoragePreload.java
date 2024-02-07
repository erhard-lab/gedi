package gedi.core.region;

import gedi.util.dynamic.DynamicObject;

public class GenomicRegionStoragePreload<T> {

	private Class<T> type;
	private T example;
	private DynamicObject metaData;

	public GenomicRegionStoragePreload(Class<T> type, T example, DynamicObject metaData) {
		this.type = type;
		this.example = example;
		this.metaData = metaData;
	}
	
	public Class<?> getType() {
		return type;
	}
	
	public T getExample() {
		return example;
	}
	
	public DynamicObject getMetaData() {
		return metaData;
	}
	
}
