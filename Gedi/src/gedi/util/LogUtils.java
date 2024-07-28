package gedi.util;
import java.io.InputStream;
import java.net.URL;
import java.util.logging.LogManager;


/**
 * Call config before any logging method is used, and then use
 * @author erhard
 *
 */
public class LogUtils extends LogManager {
	



	public enum LogMode {
		Normal,Silent,Debug
	}
	
//	public static void setFile(String path) throws SecurityException, IOException {
//		setFile(path,false);
//	}
//	
//	public static void setFile(String path, boolean append) throws SecurityException, IOException {
//		removeHandlers();
//		Logger.getLogger("global").getParent().addHandler(new ConsoleHandler(){{setLevel(Level.SEVERE);}});
//		Logger.getLogger("global").getParent().addHandler(new FileHandler(path, append));
//	}
//	
//	public static void removeHandlers() {
//		for (Handler h : Logger.getLogger("global").getParent().getHandlers())
//			Logger.getLogger("global").getParent().removeHandler(h);
//	}
	
    static LogUtils instance;
    public LogUtils() { instance = this; }
    @Override public void reset() { if (doShutdown) reset0(); }
    private void reset0() { super.reset(); }
    
    public static void shutdown() {
    	if (instance!=null) instance.reset0();
    }

    private static boolean doShutdown = true;
	
	public static void config() {
		config(LogMode.Normal,true);
	}
	
	public static void config(LogMode mode, boolean doShutdown) {
		LogUtils.doShutdown = doShutdown;
		
		String fname = System.getProperty("java.util.logging.config.file");
		if (fname == null) { // already configured, do nothing!
			try {
				String res;
				if (mode==LogMode.Silent)
					res = "/logging_silent.properties";
				else if (mode==LogMode.Debug)
					res = "/logging_debug.properties";
				else
					res = "/logging.properties";
				URL intlog = LogUtils.class.getResource(res);
				InputStream str = intlog.openStream();
				LogManager.getLogManager().readConfiguration(str);
				
				str.close();
			} catch (Exception e) {
				throw new RuntimeException("Could not read internal logging configuration!",e);
			}
			
		}
			
	}
//	
//	
//	public static void setLevel(Level level) throws SecurityException, IOException {
//		removeHandlers();
//		Logger.getLogger(Logger.GLOBAL_LOGGER_NAME).getParent().addHandler(new ConsoleHandler(){{
//			setLevel(level);
//			setFormatter(new Formatter(){
//				private final Date dat = new Date();
//				public synchronized String format(LogRecord record) {
//			        dat.setTime(record.getMillis());
//			        String source;
//			        if (record.getSourceClassName() != null) {
//			            source = record.getSourceClassName();
//			            if (record.getSourceMethodName() != null) {
//			               source += " " + record.getSourceMethodName();
//			            }
//			        } else {
//			            source = record.getLoggerName();
//			        }
//			        String message = formatMessage(record);
//			        String throwable = "";
//			        if (record.getThrown() != null) {
//			            StringWriter sw = new StringWriter();
//			            PrintWriter pw = new PrintWriter(sw);
//			            pw.println();
//			            record.getThrown().printStackTrace(pw);
//			            pw.close();
//			            throwable = sw.toString();
//			        }
//			        return String.format("%1$tY-%1$tm-%1$td %1$tH:%1$tM:%1$tS.%1$tL %4$s %5$s%6$s%n",
//			                             dat,
//			                             source,
//			                             record.getLoggerName(),
//			                             record.getLevel().getName(),
//			                             message,
//			                             throwable);
//			    }
//			});
//		}});
//	}
	
}
