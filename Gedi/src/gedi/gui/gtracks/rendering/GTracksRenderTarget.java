package gedi.gui.gtracks.rendering;

import java.awt.Color;
import java.io.IOException;

import gedi.util.io.randomaccess.BinaryWriter;

public interface GTracksRenderTarget {

	void rect(double x1, double x2, double y1, double y2, Color border, Color background);
	void text(String text, double x1, double x2, double y1, double y2, Color color, Color background);
	

	int writeRaw(BinaryWriter out, int width, int height) throws IOException;
	String getFormat();
}
