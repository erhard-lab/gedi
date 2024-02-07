package gedi.core.workspace.loader;

import java.io.IOException;
import java.nio.file.Path;

public interface WorkspaceItemLoader<T,P> {

	String[] getExtensions();

	T load(Path path) throws IOException;
	P preload(Path path) throws IOException;
	
	default PreloadInfo<T, P> getPreloadInfo(Path path) throws IOException {
		return new PreloadInfo<>(getItemClass(),preload(path));
	}

	Class<T> getItemClass();
	
	
	boolean hasOptions();
	
	void updateOptions(Path path);
	
}
