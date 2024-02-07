package gedi.util.io.text.tsv.formats;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.nio.file.Path;

public class LocationsFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<NameAnnotation>,GenomicRegionStoragePreload> {

	public static final String[] extensions = new String[]{"locations"};
	
	private boolean ignore = false;
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public MemoryIntervalTreeStorage<NameAnnotation> load(Path path)
			throws IOException {
		if (ignore)
			return new LocationFileReader(path.toString()).readIntoMemoryTakeFirst();
		else
			return new LocationFileReader(path.toString()).readIntoMemoryThrowOnNonUnique();
	}
	
	public LocationsFileLoader setIgnore(boolean ignore) {
		this.ignore = ignore;
		return this;
	}

	@Override
	public Class<MemoryIntervalTreeStorage<NameAnnotation>> getItemClass() {
		return (Class)MemoryIntervalTreeStorage.class;
	}
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		return new GenomicRegionStoragePreload(NameAnnotation.class, new NameAnnotation(""),DynamicObject.getEmpty());
	}

	@Override
	public boolean hasOptions() {
		return false;
	}

	@Override
	public void updateOptions(Path path) {
	}

}
