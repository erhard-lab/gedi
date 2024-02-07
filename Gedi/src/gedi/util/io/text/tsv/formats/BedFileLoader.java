package gedi.util.io.text.tsv.formats;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.nio.file.Path;

public class BedFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<ScoreNameAnnotation>,GenomicRegionStoragePreload> {

	public static final String[] extensions = new String[]{"bed","bed.gz"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public MemoryIntervalTreeStorage<ScoreNameAnnotation> load(Path path)
			throws IOException {
		return new ScoreNameBedFileReader(path.toString()).readIntoMemoryTakeFirst();
	}
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		return new GenomicRegionStoragePreload(ScoreNameAnnotation.class, new ScoreNameAnnotation("", 0),DynamicObject.getEmpty());
	}


	@Override
	public Class<MemoryIntervalTreeStorage<ScoreNameAnnotation>> getItemClass() {
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
