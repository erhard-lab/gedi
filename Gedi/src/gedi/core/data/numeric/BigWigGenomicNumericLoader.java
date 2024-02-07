package gedi.core.data.numeric;

import java.io.IOException;
import java.nio.file.Path;

import gedi.core.workspace.loader.WorkspaceItemLoader;

public class BigWigGenomicNumericLoader implements WorkspaceItemLoader<BigWigGenomicNumericProvider,Void> {

	public static String[] extensions = {"bigwig","bw"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public BigWigGenomicNumericProvider load(Path path) throws IOException {
		return new BigWigGenomicNumericProvider(path.toString());
	}

	@Override
	public Class<BigWigGenomicNumericProvider> getItemClass() {
		return BigWigGenomicNumericProvider.class;
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
