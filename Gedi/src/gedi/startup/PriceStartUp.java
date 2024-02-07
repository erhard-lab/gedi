package gedi.startup;

import java.lang.reflect.Modifier;

import gedi.app.Startup;
import gedi.app.classpath.ClassPath;
import gedi.app.classpath.ClassPathCache;
import gedi.core.region.GenomicRegionStorage;
import gedi.riboseq.analysis.PriceAnalysis;
import gedi.riboseq.analysis.PriceAnalysisExtensionPoint;
import gedi.riboseq.visu.PriceOrfGenomicRegionStorageDisplayTemplateGenerator;
import gedi.util.ReflectionUtils;
import gedi.util.io.text.jhp.TemplateGenerator;
import gedi.util.io.text.jhp.display.DisplayTemplateGeneratorExtensionPoint;

public class PriceStartUp implements Startup {

	@Override
	public void accept(ClassPath t) {
		DisplayTemplateGeneratorExtensionPoint.getInstance().addExtension((Class)PriceOrfGenomicRegionStorageDisplayTemplateGenerator.class, (Class)PriceOrfGenomicRegionStorageDisplayTemplateGenerator.cls);

		
		for (String cls : ClassPathCache.getInstance().getClassesOfPackage("gedi.riboseq.analysis")) {
			try {
				Class c = (Class)Class.forName("gedi.riboseq.analysis."+cls);
				if (PriceAnalysis.class.isAssignableFrom(c) 
						&& !c.isInterface() 
						&& !Modifier.isAbstract(c.getModifiers())
						&& ReflectionUtils.hasStatic(c, "name"))
					PriceAnalysisExtensionPoint.getInstance().addExtension((Class)c, (String) ReflectionUtils.getStatic(c, "name"));
			} catch (Exception e) {
				throw new RuntimeException("Could not load class "+cls,e);
			}
		}
	}

	
	
	
}
