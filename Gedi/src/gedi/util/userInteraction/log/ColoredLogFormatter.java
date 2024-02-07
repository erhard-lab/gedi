package gedi.util.userInteraction.log;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Date;
import java.util.logging.Formatter;
import java.util.logging.Level;
import java.util.logging.LogManager;
import java.util.logging.LogRecord;

public class ColoredLogFormatter extends Formatter {
    public static final String ANSI_RESET = "\u001B[0m";

    private static final String format = LogManager.getLogManager().getProperty("format");
    private final Date dat = new Date();
    
    private static final String colSevere = LogManager.getLogManager().getProperty("ansi.severe");
    private static final String colWarning = LogManager.getLogManager().getProperty("ansi.warning");
    private static final String colInfo = LogManager.getLogManager().getProperty("ansi.info");
    private static final String colOther = LogManager.getLogManager().getProperty("ansi.other");
    
    private int warnings = 0;
    private int severes = 0;
    
    
    public synchronized String format(LogRecord record) {
    	String message = formatMessage(record);
        Level lev = record.getLevel();
        
    	if (lev.intValue()==Level.OFF.intValue()) {
    		lev = Level.INFO;
    		if (warnings>0) lev = Level.WARNING;
    		if (severes>0) lev = Level.SEVERE;
    		message = "Warnings: "+warnings+", errors: "+severes+"!";
    	}
    	
        dat.setTime(record.getMillis());
        String source;
        if (record.getSourceClassName() != null) {
            source = record.getSourceClassName();
            if (record.getSourceMethodName() != null) {
               source += " " + record.getSourceMethodName();
            }
        } else {
            source = record.getLoggerName();
        }
        String throwable = "";
        if (record.getThrown() != null) {
            StringWriter sw = new StringWriter();
            PrintWriter pw = new PrintWriter(sw);
            pw.println();
            record.getThrown().printStackTrace(pw);
            pw.close();
            throwable = sw.toString();
        }
        String col;
        if (lev.intValue()==Level.SEVERE.intValue()) col = colSevere;
        else if (lev.intValue()==Level.WARNING.intValue()) col = colWarning;
        else if (lev.intValue()==Level.INFO.intValue()) col = colInfo;
        else col = colOther;
        if (col==null) col=ANSI_RESET;
        
        if (lev.intValue()==Level.SEVERE.intValue()) severes++;
        if (lev.intValue()==Level.WARNING.intValue()) warnings++;
        
        return col+String.format(format,
                             dat,
                             source,
                             record.getLoggerName(),
                             lev.getLocalizedName()+ANSI_RESET,
                             message,
                             throwable);
    }
}