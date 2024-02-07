package gedi.util.program;

import java.io.File;
import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.function.Supplier;
import java.util.logging.Logger;

import gedi.util.StringUtils;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.userInteraction.progress.Progress;

public class GediProgramContext {

	
	private Logger log;
	private Supplier<Progress> progress;
	private boolean dry;
	private GediParameterSet params;
	private File runFile;
	
	private Date start;
	private HashMap<String,Date> starts = new HashMap<String, Date>();   
			
	
	public GediProgramContext(Logger log, Supplier<Progress> progress, File runFile, boolean dry, GediParameterSet params) {
		this.log = log;
		this.progress = progress;
		this.dry = dry;
		this.runFile = runFile;
		this.params = params;
	}

	public <T extends GediParameterSet> T getParams() {
		return (T) params;
	}
	
	public Logger getLog() {
		return log;
	}
	
	public File getRunFile() {
		return runFile;
	}
	
	public void startMessage(String name) throws IOException {
		String msg = String.format("Starting %s ...",name);
		log.info(msg);
		starts.put(name,new Date());
		if (runFile!=null) 
			writeIntoRunFile(msg);
	}
	
	public void finishMessage(String name) throws IOException {
		String span = starts.containsKey(name)?StringUtils.getHumanReadableTimespan(new Date().getTime()-starts.get(name).getTime()):"unknown";
		String msg = String.format("Finished %s (time: %s)",name,span);
		log.info(msg);
		if (runFile!=null) 
			writeIntoRunFile(msg);
	}
	
	public void writeIntoRunFile(String status) throws IOException {
		if (runFile==null) return;
		
		if (start==null) start=new Date();
		Date now = new Date();
		
		synchronized (runFile) {
			LineWriter out = new LineOrientedFile(runFile.getPath()).append();
	        String msg = String.format("%1$s\t%2$tY-%2$tm-%2$td %2$tH:%2$tM:%2$tS.%2$tL %3$s%n",
	        		StringUtils.getHumanReadableTimespan(now.getTime()-start.getTime()),
                    now,
                    status);
			out.write(msg);
			out.close();
		}
	}
		
	public void writeIntoRunFile(String format, Object...args) throws IOException {
		writeIntoRunFile(String.format(format, args));
	}

	public Progress getProgress() {
		return progress.get();
	}
	
	public boolean isDryRun() {
		return dry;
	}
	
	public void logf(String format, Object...args) {
		getLog().info(String.format(format, args));
	}
	
	public void errorf(String format, Object...args) {
		getLog().severe(String.format(format, args));
	}
	
	public void warningf(String format, Object...args) {
		getLog().warning(String.format(format, args));
	}
	
	
}
