package gedi.gui.gtracks.rendering.target;

import java.awt.Dimension;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.freehep.graphicsio.svg.SVGGraphics2D;

import gedi.util.StringUtils;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class SvgRenderTarget extends Graphics2DRenderTarget<SVGGraphics2D> {

	
	private ByteArrayOutputStream buff;


	public SvgRenderTarget() {
		buff = new ByteArrayOutputStream(1024*64);
		g2 = new SVGGraphics2D(buff, new Dimension(1,1));
		g2.startExport();
	}
	
	@Override
	public String getFormat() {
		return "svg";
	} 
	

	private String data;
	public int writeRaw(BinaryWriter out, int width, int height) throws IOException {
		if (data==null)  {
			g2.endExport();
			
			// clean up and insert width/height
			String src = new String(buff.toByteArray());
			src = src.substring(src.indexOf("</desc>"));
			
			src = String.format(
					"<svg version=\"1.1\" baseProfile=\"full\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
			+ "xmlns:ev=\"http://www.w3.org/2001/xml-events\" x=\"0px\" y=\"0px\" width=\"%dpx\" height=\"%dpx\" viewBox=\"0 0 %d %d\"><desc>Creator: Gedi",
			width,height,width,height)+src;
			data = src;
		}
		out.putAsciiChars(data);
		return data.length();
	}
	
}
