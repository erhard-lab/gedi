package gedi.util.io.text.tsv.formats;

import java.io.IOException;
import java.nio.file.Path;

import gedi.core.data.annotation.NarrowPeakAnnotation;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

public class NarrowPeakFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<NarrowPeakAnnotation>,GenomicRegionStoragePreload> {

	public static final String[] extensions = new String[]{"narrowPeak","narrowPeak.gz"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public MemoryIntervalTreeStorage<NarrowPeakAnnotation> load(Path path)
			throws IOException {
		return new NarrowPeakFileReader(path.toString()).readIntoMemoryTakeFirst();
	}
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		return new GenomicRegionStoragePreload(NarrowPeakAnnotation.class, new NarrowPeakAnnotation(),DynamicObject.getEmpty());
	}


	@Override
	public Class<MemoryIntervalTreeStorage<NarrowPeakAnnotation>> getItemClass() {
		return (Class)MemoryIntervalTreeStorage.class;
	}

	@Override
	public boolean hasOptions() {
		return false;
	}

	@Override
	public void updateOptions(Path path) {
	}

}
