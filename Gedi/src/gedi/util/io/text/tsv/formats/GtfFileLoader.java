package gedi.util.io.text.tsv.formats;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.nio.file.Path;

public class GtfFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<Transcript>,GenomicRegionStoragePreload> {

	public static final String[] extensions = new String[]{"gtf"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public MemoryIntervalTreeStorage<Transcript> load(Path path)
			throws IOException {
		return new GtfFileReader(path.toString(), "exon", "CDS").readIntoMemoryThrowOnNonUnique();
	}
	
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		return new GenomicRegionStoragePreload(Transcript.class, new Transcript(),DynamicObject.getEmpty());
	}


	@Override
	public Class<MemoryIntervalTreeStorage<Transcript>> getItemClass() {
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
