package gedi.util.userInteraction.progress;

import java.io.PrintStream;
import java.util.Locale;

import gedi.util.StringUtils;
import gedi.util.functions.EI;

public class ConsoleProgress extends AbstractProgress {

	private int pad = 0;
	private int maxDescrLen = 80;
	public PrintStream out;
	
	
	public ConsoleProgress() {
		this(System.out);
	}

	public ConsoleProgress(PrintStream out) {
		this(out,null);
	}

	public ConsoleProgress(ProgressManager man) {
		this(System.out,man);
	}

	public ConsoleProgress(PrintStream out, ProgressManager man) {
		this.out = out;
		if (man!=null)
			setManager(man);
	}
	
	public ConsoleProgress setTheme(int theme) {
		this.theme = theme&viewcss.length;
		return this;
	}

	public void setMaxDescriptionLength(int maxDescrLen) {
		this.maxDescrLen = maxDescrLen;
	}

	private int theme = 4;
	private int view;
	private String[][] viewcss = {
				{"|","/","-","\\"},
				EI.substrings("▁▃▄▅▆▇█▇▆▅▄▃",1).toArray(String.class),
				{"▹▹▹▹▹","▸▹▹▹▹","▹▸▹▹▹","▹▹▸▹▹","▹▹▹▸▹","▹▹▹▹▸"},
				EI.substrings("▘▀▜█▟▄▖▌▛█▜▐▗▄▙█▛▀▝▐▟█▙▌",1).toArray(String.class),
				{"▀▘","▝▀"," ▜"," ▟","▗▄","▄▖","▙ ","▛ "},
				{"▀  "," ▀ ","  ▀","  ▐","  ▄"," ▄ ","▄  ","▌  ","▀  ","▀▀ ","▀▀▀","▀▀▜","▀▀█","▀██","███","▟██","▄██","▄▄█","▄▄▟","▄▄▄","▄▄ ","▄  ","   "}
	};
			
	private static boolean first = true;

	@Override
	public void firstView(int total) {
		if (total>1)
			out.println();
		if (total==1) {
			out.print("\33[?25l");
			if (first) {
				Runtime.getRuntime().addShutdownHook(new Thread(()->out.print("\33[?25h")));
				first = false;
			}
		}
	}

	@Override
	public void updateView(int index, int total) {
		CharSequence d = getDescription();
		if (d==null) d="";
		if (d.length()>maxDescrLen) d = d.subSequence(0, maxDescrLen-3).toString()+"...";
		if (d.length()>pad) pad = d.length();
		
		String description = StringUtils.padRight(d.toString(), pad,' ');
		
		view++;
		
		if (index==0 && total>1) {
			out.print("\33["+(total-1)+"A");
		}
		
		String viewc = viewcss[theme][view%viewcss[theme].length];
		
		if (progress>0) {
			if (isGoalKnown())
				out.printf(Locale.US,"\33[2K\u001B[44m\u001B[1m%s\u001B[0m %s %d/%d (%.2f%%, %.1f/sec) Estimated time: %s\r",viewc,description,progress,count,(progress*100.0)/count,getPerSecond(),getEstimatedTime());
			else 
				out.printf(Locale.US,"\33[2K\u001B[44m\u001B[1m%s\u001B[0m %s %d (%.1f/sec)\r",viewc,description,progress,getPerSecond());
		}
		else {
			out.printf(Locale.US,"\33[2K\u001B[44m\u001B[1m%s\u001B[0m %s\r",viewc,description);
		}
		
		if (index+1<total)
			out.print("\n");
	}
	
	
	@Override
	public void lastView(int total) {
		if (total>1)
			out.print("\33["+(total-1)+"A");
		out.printf(Locale.US,"\r\33[2KProcessed %d elements in %s (Throughput: %.1f/sec)\n",progress,StringUtils.getHumanReadableTimespan(getTotalTime()),getPerSecond());
		if (total>1)
			out.print("\33["+(total-1)+"B");
		if (total==1)
			out.print("\33[?25h");
	}

}
