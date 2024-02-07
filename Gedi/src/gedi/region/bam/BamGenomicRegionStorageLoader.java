package gedi.region.bam;

import gedi.core.data.HasConditions;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.util.dynamic.DynamicObject;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.nashorn.JS;

import java.io.IOException;
import java.nio.file.Path;

import javax.script.ScriptException;

public class BamGenomicRegionStorageLoader implements WorkspaceItemLoader<BamGenomicRegionStorage,GenomicRegionStoragePreload> {

	public static String[] extensions = {"bam","bamlist"};
	
	@Override
	public String[] getExtensions() {
		return extensions;
	}

	@Override
	public BamGenomicRegionStorage load(Path path) throws IOException {
		String p = path.toString();
		if (p.endsWith(".bam"))
			return new BamGenomicRegionStorage(p).detectVariations();
		else if (p.endsWith(".bamlist")) {
			BamGenomicRegionStorage re = new BamGenomicRegionStorage(new LineOrientedFile(p).readAllLines("#")).detectVariations();
			String src = new LineOrientedFile(p).lineIterator().filterString(s->s.startsWith("#!")).substring(2).concat("\n");
			if (src.length()>0)
				try {
					new JS().putVariable("bam", re).execSource(src);
				} catch (ScriptException e) {
					throw new IOException("Cannot execute commands in bamlist!", e);
				}
			return re;
		}
		throw new RuntimeException("Unknown file type: "+p);
	}

	@Override
	public Class<BamGenomicRegionStorage> getItemClass() {
		return BamGenomicRegionStorage.class;
	}

	@Override
	public GenomicRegionStoragePreload preload(Path path) throws IOException {
		BamGenomicRegionStorage bam = load(path);
		Class<?> cls = bam.getType();
		AlignedReadsData rec = bam.getRandomRecord();
		DynamicObject meta = bam.getMetaData();
		bam.close();
		return new GenomicRegionStoragePreload(AlignedReadsData.class, rec, meta);
	}
	
	@Override
	public boolean hasOptions() {
		return false;
	}

	@Override
	public void updateOptions(Path path) {
		
	}

}
