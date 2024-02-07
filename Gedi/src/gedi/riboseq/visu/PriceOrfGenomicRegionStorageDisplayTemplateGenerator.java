package gedi.riboseq.visu;

import java.io.File;

import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.NameProvider;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.workspace.loader.PreloadInfo;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.startup.CoreStartup;
import gedi.util.functions.EI;
import gedi.util.io.text.jhp.TemplateEngine;
import gedi.util.io.text.jhp.TemplateGenerator;
import gedi.util.io.text.jhp.display.AutoTemplateGenerator;
import gedi.util.io.text.tsv.formats.BedFileLoader;

public class PriceOrfGenomicRegionStorageDisplayTemplateGenerator implements AutoTemplateGenerator<PreloadInfo<?, GenomicRegionStoragePreload>,PriceOrfGenomicRegionStorageDisplayTemplateGenerator> {

	
	/**
	 * For auto inclusion in {@link CoreStartup}
	 */
	public static final Class<?> cls = GenomicRegionStorage.class;
	
	
	private String path;
	
	public PriceOrfGenomicRegionStorageDisplayTemplateGenerator(String path) {
		this.path = path;
	}

	@Override
	public void accept(TemplateEngine t, PriceOrfGenomicRegionStorageDisplayTemplateGenerator[] g) {
		Runnable res = t.save("paths","id");
		t.set("paths", EI.wrap(g).map(te->te.path).toArray(String.class));
		t.set("id", path);
		t.template("priceorfgenomicregion");
		res.run();
	}

	@Override
	public double applyAsDouble(PreloadInfo<?, GenomicRegionStoragePreload> pre) {
		GenomicRegionStoragePreload value = pre.getInfo();
		if (PriceOrf.class.isAssignableFrom(value.getType())) return 10;
		return -1;
	}

}
