package gedi.core.data.numeric.diskrmq;

import java.io.IOException;
import java.nio.file.Path;

import gedi.core.workspace.loader.WorkspaceItemLoader;

public class DiskGenomicNumericLoader implements WorkspaceItemLoader<DiskGenomicNumericProvider,Void> {

	public static String[] extensions = {"rmq"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public DiskGenomicNumericProvider load(Path path) throws IOException {
		return new DiskGenomicNumericProvider(path.toString());
	}

	@Override
	public Class<DiskGenomicNumericProvider> getItemClass() {
		return DiskGenomicNumericProvider.class;
	}
	
	@Override
	public Void preload(Path path) throws IOException {
		return null;
	}

	@Override
	public boolean hasOptions() {
		return false;
	}

	@Override
	public void updateOptions(Path path) {
		
	}

}
