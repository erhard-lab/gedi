package gedi.gui.gtracks.rendering;

@FunctionalInterface
public interface GTracksRenderer {

	
	GTracksRenderRequest render(GTracksRenderTarget target, GTracksRenderContext context);
	
}
