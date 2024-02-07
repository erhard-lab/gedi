package gedi.util.io.text.tsv.formats;

import gedi.core.data.annotation.Gff3Element;
import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.nio.file.Path;

public class GffFileLoader implements WorkspaceItemLoader<MemoryIntervalTreeStorage<NameAnnotation>,GenomicRegionStoragePreload> {

	public static final String[] extensions = new String[]{"gff","gff3"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public MemoryIntervalTreeStorage<NameAnnotation> load(Path path)
			throws IOException {
		return new GffFileReader<>(path.toString(), gff->new NameAnnotation(decideName(gff)), NameAnnotation.class).readIntoMemoryThrowOnNonUnique();
	}
	
	private static String decideName(Gff3Element gff) {
		if (gff.getAttributeNames().contains("Name"))
			return (String) gff.getAttribute("Name");
		
		if (gff.getAttributeNames().contains("name"))
			return (String) gff.getAttribute("name");
		
		if (gff.getAttributeNames().contains("ID"))
			return (String) gff.getAttribute("ID");
		
		if (gff.getAttributeNames().contains("id"))
			return (String) gff.getAttribute("id");
		
		return gff.toString();
	}
	
	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		return new GenomicRegionStoragePreload(NameAnnotation.class, new NameAnnotation(""),DynamicObject.getEmpty());
	}


	@Override
	public Class<MemoryIntervalTreeStorage<NameAnnotation>> getItemClass() {
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
