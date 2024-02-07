package gedi.gui.gtracks.rendering.target;

import java.awt.Color;
import java.io.IOException;

import gedi.gui.gtracks.rendering.GTracksRenderTarget;
import gedi.util.io.randomaccess.BinaryWriter;

public class CompositeRenderTarget<A extends GTracksRenderTarget,B extends GTracksRenderTarget> implements GTracksRenderTarget{

	private A a;
	private B b;
	
	public CompositeRenderTarget(A a, B b) {
		this.a = a;
		this.b = b;
	}

	@Override
	public void rect(double x1, double x2, double y1, double y2, Color border, Color background) {
		a.rect(x1, x2, y1, y2, border, background);
		b.rect(x1, x2, y1, y2, border, background);
	}

	@Override
	public void text(String text, double x1, double x2, double y1, double y2, Color color, Color background) {
		a.text(text, x1, x2, y1, y2, color, background);
		b.text(text, x1, x2, y1, y2, color, background);
	}
	
	public A a() { return a; }
	public B b() { return b; }


	@Override
	public int writeRaw(BinaryWriter out, int width, int height) throws IOException {
		return a.writeRaw(out, width, height);
	}
	
	@Override
	public String getFormat() {
		return a.getFormat();
	}
}
