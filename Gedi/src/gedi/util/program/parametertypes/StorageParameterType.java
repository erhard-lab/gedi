package gedi.util.program.parametertypes;


import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.workspace.loader.WorkspaceItemLoaderExtensionPoint;

public class StorageParameterType<T> implements GediParameterType<GenomicRegionStorage<T>> {

	@Override
	public GenomicRegionStorage<T> parse(String s) {
		return GenomicRegionStorage.load(s);
	}

	@Override
	public Class<GenomicRegionStorage<T>> getType() {
		return (Class)GenomicRegionStorage.class;
	}
	

}
