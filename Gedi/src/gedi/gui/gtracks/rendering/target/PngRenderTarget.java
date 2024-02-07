package gedi.gui.gtracks.rendering.target;

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import gedi.util.PaintUtils;
import gedi.util.io.randomaccess.BinaryWriter;

public class PngRenderTarget extends Graphics2DRenderTarget<Graphics2D> {

	private static final double fac = 1.5;
	
	private BufferedImage img;
	
	public PngRenderTarget() {
		img = new BufferedImage(1920, 400, BufferedImage.TYPE_INT_ARGB);
		g2 = img.createGraphics();
	}
	
	@Override
	public String getFormat() {
		return "png";
	}
	
	@Override
	protected void checkBoundindBox(double width, double height) {
		int nw = img.getWidth();
		while (nw<width) nw = (int) Math.ceil(nw*fac);
		int nh = img.getHeight();
		while (nh<height) nh = (int) Math.ceil(nh*fac);
		if (nh!=img.getHeight() || nw!=img.getWidth()) {
			g2.dispose();
			img = PaintUtils.resize(img,nw,nh);
			g2 = img.createGraphics();
		}
	}

	public BufferedImage getImage() {
		return img;
	}

	private byte[] bbuf; 

	@Override
	public int writeRaw(BinaryWriter out, int width, int height) throws IOException {
		if (bbuf==null) {
			ImageIO.write(img, "png", new File("before.png"));
			if (height!=img.getHeight() || width!=img.getWidth())
				img = PaintUtils.resize(img,width,height);
			ImageIO.write(img, "png", new File("after.png"));
			ByteArrayOutputStream buff = new ByteArrayOutputStream();
			ImageIO.write(img, "png", buff);
			bbuf = buff.toByteArray();
			
			g2.dispose();
			g2 = null;
		}
		out.put(bbuf, 0, bbuf.length);
		return bbuf.length;
	}
	
}
