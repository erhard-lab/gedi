package gedi.util.io.text.jhp.display;

import java.io.File;

import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.workspace.loader.PreloadInfo;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.startup.CoreStartup;
import gedi.util.functions.EI;
import gedi.util.io.text.jhp.TemplateEngine;
import gedi.util.io.text.jhp.TemplateGenerator;
import gedi.util.io.text.tsv.formats.BedFileLoader;

public class ScoreNameGenomicRegionStorageDisplayTemplateGenerator implements AutoTemplateGenerator<PreloadInfo<?, GenomicRegionStoragePreload>,ScoreNameGenomicRegionStorageDisplayTemplateGenerator> {

	
	/**
	 * For auto inclusion in {@link CoreStartup}
	 */
	public static final Class<?> cls = GenomicRegionStorage.class;
	
	
	private String path;
	
	public ScoreNameGenomicRegionStorageDisplayTemplateGenerator(String path) {
		this.path = path;
	}

	@Override
	public void accept(TemplateEngine t, ScoreNameGenomicRegionStorageDisplayTemplateGenerator[] g) {
		Runnable res = t.save("paths","id");
		t.set("paths", EI.wrap(g).map(te->te.path).toArray(String.class));
		t.set("id", path);
		t.template("scorenamegenomicregion");
		res.run();
	}

	@Override
	public double applyAsDouble(PreloadInfo<?, GenomicRegionStoragePreload> pre) {
		GenomicRegionStoragePreload value = pre.getInfo();
		if (value.getType()==ScoreNameAnnotation.class) return 10;
		if (ScoreNameAnnotation.class.isAssignableFrom(value.getType())) return 1;
		return -1;
	}

}
