package gedi.gui.gtracks.rendering.target;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import gedi.gui.gtracks.rendering.GTracksRenderTarget;
import gedi.util.PaintUtils;

public abstract class Graphics2DRenderTarget<G extends Graphics2D> implements GTracksRenderTarget {

	
	protected G g2;

	protected void checkBoundindBox(double width, double height) {
	}

	@Override
	public void rect(double x1, double x2, double y1, double y2, Color border, Color background) {
		checkBoundindBox(x2, y2);
		if (background!=null) {
			g2.setColor(background);
			if (border==null)
				g2.draw(new Rectangle2D.Double(x1, y1, x2-x1, y2-y1));
			g2.fill(new Rectangle2D.Double(x1, y1, x2-x1, y2-y1));
		}
		if (border!=null) {
			g2.setColor(border);
			g2.draw(new Rectangle2D.Double(x1, y1, x2-x1, y2-y1));
		}
	}

	@Override
	public void text(String text, double x1, double x2, double y1, double y2, Color color, Color background) {
		checkBoundindBox(x2, y2);
		if (background!=null)
			rect(x1, x2, y1, y2, null, background);
		g2.setColor(color);
		PaintUtils.paintString(text, g2, new Rectangle2D.Double(x1, y1, x2-x1, y2-y1), 0, 0);
	}
	
}
