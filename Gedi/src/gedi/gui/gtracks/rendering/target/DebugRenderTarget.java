package gedi.gui.gtracks.rendering.target;

import java.awt.Color;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;

import gedi.gui.gtracks.rendering.GTracksRenderTarget;
import gedi.util.io.randomaccess.BinaryWriter;

public class DebugRenderTarget implements GTracksRenderTarget {

	private static Logger log = Logger.getLogger(DebugRenderTarget.class.getName());
	
	private ArrayList<String> commands = new ArrayList<>(); 
	
	@Override
	public String getFormat() {
		return "debug";
	}
	
	@Override
	public void rect(double x1, double x2, double y1, double y2, Color border, Color background) {
		String s = String.format("rect(%.2f,%.2f,%.2f,%.2f,%s,%s)",x1,x2,y1,y2,border,background);
		log.info(s);
		commands.add(s);
	}

	@Override
	public void text(String text, double x1, double x2, double y1, double y2, Color color, Color background) {
		String s = String.format("text(%s,%.2f,%.2f,%.2f,%.2f,%s,%s)",text,x1,x2,y1,y2,color,background);	
		log.info(s);
		commands.add(s);	
	}

	public ArrayList<String> getCommands() {
		return commands;
	}
	
	public String[] getAndClearCommands() {
		String[] re = commands.toArray(new String[0]);
		commands.clear();
		return re;
	}

	@Override
	public int writeRaw(BinaryWriter out, int width, int height) throws IOException {
		int re = 0;
		for (String c : commands) {
			re+=c.length()+1;
			out.putAsciiChars(c);
			out.putAsciiChar('\n');
		}
		return re;
	}
	
}
