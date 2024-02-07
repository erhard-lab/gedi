package gedi.util.io.text.jhp;

import java.io.File;
import java.io.IOException;

import gedi.app.classpath.ClassPathCache;
import gedi.util.StringUtils;
import gedi.util.io.text.LineIterator;

public class TemplateExtensions {

	private Jhp jhp;
	

	public TemplateExtensions(Jhp jhp) {
		this.jhp = jhp;
	}

	@SuppressWarnings("resource")
	public void include(String id) throws IOException {
		String inc = new LineIterator(TemplateExtensions.class.getResource(id).openStream()).concat("\n");
		jhp.getJs().getStdout().write(inc);
	}
	
	
	public void includeExtensions(String point) throws IOException {
		File f = new File(point);
		for (String s : ClassPathCache.getInstance().getResourcesOfPath(StringUtils.removeHeader(f.getParent(),"/")))
			if (s.startsWith(f.getName()+"-"))
				include(f.getParent()+"/"+s);
	}
	
}
