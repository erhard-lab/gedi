package gedi.centeredDiskIntervalTree;

import java.io.IOException;
import java.nio.file.Path;

import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

public class CenteredDiskIntervalTreeStorageLoader implements WorkspaceItemLoader<CenteredDiskIntervalTreeStorage,GenomicRegionStoragePreload> {

	private static String[] extensions = {"cit"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public CenteredDiskIntervalTreeStorage<?> load(Path path) throws IOException {
		return new CenteredDiskIntervalTreeStorage(path.toString());
	}

	@Override
	public Class<CenteredDiskIntervalTreeStorage> getItemClass() {
		return CenteredDiskIntervalTreeStorage.class;
	}
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		CenteredDiskIntervalTreeStorage cit = new CenteredDiskIntervalTreeStorage(path.toString());
		Class<?> cls = cit.getType();
		Object rec = cit.getRandomRecord();
		DynamicObject meta = cit.getMetaData();
		cit.close();
		return new GenomicRegionStoragePreload(cls, rec, meta);
	}

	@Override
	public boolean hasOptions() {
		return false;
	}

	@Override
	public void updateOptions(Path path) {
		
	}

}
