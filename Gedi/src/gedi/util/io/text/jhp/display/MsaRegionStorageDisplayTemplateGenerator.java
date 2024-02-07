package gedi.util.io.text.jhp.display;

import java.io.File;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.NameProvider;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.workspace.loader.PreloadInfo;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.startup.CoreStartup;
import gedi.util.algorithm.string.alignment.multiple.MsaBlock;
import gedi.util.functions.EI;
import gedi.util.io.text.jhp.TemplateEngine;
import gedi.util.io.text.jhp.TemplateGenerator;
import gedi.util.io.text.tsv.formats.BedFileLoader;

public class MsaRegionStorageDisplayTemplateGenerator implements AutoTemplateGenerator<PreloadInfo<?, GenomicRegionStoragePreload>,MsaRegionStorageDisplayTemplateGenerator> {

	
	/**
	 * For auto inclusion in {@link CoreStartup}
	 */
	public static final Class<?> cls = GenomicRegionStorage.class;
	
	
	private String path;
	
	public MsaRegionStorageDisplayTemplateGenerator(String path) {
		this.path = path;
	}

	@Override
	public void accept(TemplateEngine t, MsaRegionStorageDisplayTemplateGenerator[] g) {
		Runnable res = t.save("path","id");
		if (g.length>1) throw new RuntimeException("Only a single MSA allowed!");
		t.set("path", g[0].path);
		t.set("id", path);
		t.template("msagenomicregion");
		res.run();
	}

	@Override
	public double applyAsDouble(PreloadInfo<?, GenomicRegionStoragePreload> pre) {
		GenomicRegionStoragePreload value = pre.getInfo();
		if (value.getType()==MsaBlock.class) return 10;
		if (MsaBlock.class.isAssignableFrom(value.getType())) return 1;
		return -1;
	}

}
