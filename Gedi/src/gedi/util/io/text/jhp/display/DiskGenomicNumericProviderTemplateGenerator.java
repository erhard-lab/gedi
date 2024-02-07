package gedi.util.io.text.jhp.display;

import gedi.core.data.numeric.diskrmq.DiskGenomicNumericProvider;
import gedi.core.workspace.loader.PreloadInfo;
import gedi.startup.CoreStartup;
import gedi.util.functions.EI;
import gedi.util.io.text.jhp.TemplateEngine;

public class DiskGenomicNumericProviderTemplateGenerator implements AutoTemplateGenerator<PreloadInfo<?, Void>,DiskGenomicNumericProviderTemplateGenerator> {

	
	/**
	 * For auto inclusion in {@link CoreStartup}
	 */
	public static final Class<?> cls = DiskGenomicNumericProvider.class;
	
	
	private String path;
	
	public DiskGenomicNumericProviderTemplateGenerator(String path) {
		this.path = path;
	}

	@Override
	public void accept(TemplateEngine t, DiskGenomicNumericProviderTemplateGenerator[] g) {
		Runnable res = t.save("paths","id");
		t.set("paths", EI.wrap(g).map(te->te.path).toArray(String.class));
		t.set("id", path);
		t.template("diskgenomicnumericprovider");
		res.run();
	}

	@Override
	public double applyAsDouble(PreloadInfo<?, Void> pre) {
		return 1;
	}

}
