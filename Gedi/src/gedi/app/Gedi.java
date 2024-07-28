package gedi.app;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.Locale;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;

import gedi.app.classpath.ClassPath;
import gedi.app.classpath.ClassPathCache;
import gedi.app.classpath.JARClassPath;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.core.workspace.loader.WorkspaceItemLoaderExtensionPoint;
import gedi.util.LogUtils;
import gedi.util.LogUtils.LogMode;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

public class Gedi {
	
	static {
        // must be called before any Logger method is used.
        System.setProperty("java.util.logging.manager", LogUtils.class.getName());
    }
	
	private static boolean started = false;
	
	public static void startup() {
		startup(false);
	}
	public synchronized static void startup(boolean discoverClasses) {
		startup(discoverClasses,isDebug()?LogMode.Debug:LogMode.Normal, null);
	}
	private static boolean isDebug() {
		return java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0;
	}
	public synchronized static void startup(boolean discoverClasses, LogMode mode, String app) {
		if (!started) {
			Locale.setDefault(Locale.US);
			long start = System.currentTimeMillis();
			Config.getInstance();
			LogUtils.config(mode,false);
			
			Logger log = getLog();
			log.info("OS: "+System.getProperty("os.name")+" "+System.getProperty("os.version")+" "+System.getProperty("os.arch"));
			log.info("Java: "+System.getProperty("java.vm.name")+" "+System.getProperty("java.runtime.version"));
			log.info("Gedi version "+version()+" ("+develOption()+") startup");
			if (isDebug()) 
				log.info("Debug mode");
			
			log.info("Command: "+getStartupCommand());
			
			if (app!=null)
				log.info(app+" version "+version(app));
			
			Runtime.getRuntime().addShutdownHook(new Thread(new ShutdownRunnable(log,start)));
			
			if (discoverClasses)
				ClassPathCache.getInstance().discover();
			ClassPathCache.getInstance().startup();
		}
		else if (discoverClasses)
			ClassPathCache.getInstance().discover();
		
		started = true;
	}
	
	private static class ShutdownRunnable implements Runnable {

		private Logger log;
		private long start;
		
		public ShutdownRunnable(Logger log, long start) {
			this.log = log;
			this.start = start;
		}

		@Override
		public void run() {
//			System.out.println(Thread.getAllStackTraces().keySet());
			log.info("Finished: "+getStartupCommand());
			log.log(Level.OFF,"");
			log.info("Took "+StringUtils.getHumanReadableTimespan(System.currentTimeMillis()-start));
			LogUtils.shutdown();
		}
		
	}
	
	public static String getStartupCommand() {
		String cmd = System.getProperty("sun.java.command");
		if (cmd==null) return "<EMBEDDED>";
		if (cmd.startsWith("executables."))
			return "gedi -e "+cmd.substring("executables.".length());
		return "gedi "+cmd;
	}
	
	private static Properties versions;
	private synchronized static Properties getVersions() throws IOException {
		if (versions==null) {
			versions = new Properties();
			for (ClassPath cp : ClassPathCache.getInstance().getClasspath()) {
				for (String res : cp.listResources("resources"))
					if (res.endsWith(".version"))
						versions.load(cp.getResourceAsStream("/resources/"+res));
			}
		}
		return versions;
	}
	
	public static String develOption() {
		return ClassPathCache.getInstance().getClassPathOfClass(Gedi.class) instanceof JARClassPath?"JAR":"devel";
	}
	
	public static String version(String app) {
		try {
			return getVersions().getProperty(app);
		} catch (IOException e) {
			throw new RuntimeException("Could not read version!",e);
		}
	}
	
	public static String version() {
		return version("Gedi");
	}
	
	public static List<String> apps() {
		try {
			return EI.wrap(getVersions().keySet()).cast(String.class).list();
		} catch (IOException e) {
			throw new RuntimeException("Could not read versions!",e);
		}
	}
	
	
	public static <T> T load(String path) throws IOException {
		WorkspaceItemLoader<T,?> loader = WorkspaceItemLoaderExtensionPoint.getInstance().get(Paths.get(path));
		if (loader==null) throw new RuntimeException("No loader for "+path);
		return loader.load(Paths.get(path));
	}
	
	public static Logger getLog() {
		return Logger.getLogger("GEDI");
	}
	
	public static void logf(String format, Object...args) {
		getLog().info(String.format(format, args));
	}
	
	public static void errorf(String format, Object...args) {
		getLog().severe(String.format(format, args));
	}
	
	public static void warningf(String format, Object...args) {
		getLog().warning(String.format(format, args));
	}
	
}
