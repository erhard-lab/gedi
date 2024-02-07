package gedi.startup;

import java.lang.reflect.Modifier;

import gedi.app.Startup;
import gedi.app.classpath.ClassPath;
import gedi.app.classpath.ClassPathCache;
import gedi.core.data.numeric.BigWigGenomicNumericLoader;
import gedi.core.data.numeric.diskrmq.DiskGenomicNumericLoader;
import gedi.core.region.GenomicRegionStorageCapabilities;
import gedi.core.region.GenomicRegionStorageExtensionPoint;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoaderExtensionPoint;
import gedi.util.ReflectionUtils;
import gedi.util.io.text.jhp.TemplateGenerator;
import gedi.util.io.text.jhp.display.DisplayTemplateGeneratorExtensionPoint;
import gedi.util.io.text.tsv.formats.BedFileLoader;
import gedi.util.io.text.tsv.formats.GediViewFileLoader;
import gedi.util.io.text.tsv.formats.GffFileLoader;
import gedi.util.io.text.tsv.formats.LocationsFileLoader;
import gedi.util.io.text.tsv.formats.NarrowPeakFileLoader;
import gedi.util.job.pipeline.ClusterPipelineRunner;
import gedi.util.job.pipeline.ParallelPipelineRunner;
import gedi.util.job.pipeline.PipelineRunnerExtensionPoint;
import gedi.util.job.pipeline.SerialPipelineRunner;
import gedi.util.oml.OmlLoader;

public class CoreStartup implements Startup {

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public void accept(ClassPath t) {
		
		for (String ext : LocationsFileLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(LocationsFileLoader.class,ext);
		
		for (String ext : BedFileLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(BedFileLoader.class,ext);

		for (String ext : NarrowPeakFileLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(NarrowPeakFileLoader.class,ext);

		for (String ext : GffFileLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(GffFileLoader.class,ext);
		
		for (String ext : GediViewFileLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(GediViewFileLoader.class,ext);
		
		for (String ext : DiskGenomicNumericLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(DiskGenomicNumericLoader.class,ext);
		
		for (String ext : BigWigGenomicNumericLoader.extensions)
			WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(BigWigGenomicNumericLoader.class,ext);
		
		
		GenomicRegionStorageExtensionPoint.getInstance().addExtension(MemoryIntervalTreeStorage.class,
				GenomicRegionStorageCapabilities.Memory,
				GenomicRegionStorageCapabilities.Add,
				GenomicRegionStorageCapabilities.Fill
				);
		
		PipelineRunnerExtensionPoint.getInstance().addExtension(SerialPipelineRunner.class, SerialPipelineRunner.name);
		PipelineRunnerExtensionPoint.getInstance().addExtension(ParallelPipelineRunner.class, ParallelPipelineRunner.name);
		PipelineRunnerExtensionPoint.getInstance().addExtension(ClusterPipelineRunner.class, ClusterPipelineRunner.name);
		
		
		
		WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(OmlLoader.class,"oml");
		WorkspaceItemLoaderExtensionPoint.getInstance().addExtension(OmlLoader.class,"oml.jhp");
		
		
		for (String cls : ClassPathCache.getInstance().getClassesOfPackage("gedi.util.io.text.jhp.display")) {
			try {
				Class c = (Class)Class.forName("gedi.util.io.text.jhp.display."+cls);
				if (TemplateGenerator.class.isAssignableFrom(c) 
						&& !c.isInterface() 
						&& !Modifier.isAbstract(c.getModifiers())
						&& ReflectionUtils.hasStatic(c, "cls"))
					DisplayTemplateGeneratorExtensionPoint.getInstance().addExtension((Class)c, (Class) ReflectionUtils.getStatic(c, "cls"));
			} catch (Exception e) {
				throw new RuntimeException("Could not load class "+cls,e);
			}
		}
	}
	
	
}
