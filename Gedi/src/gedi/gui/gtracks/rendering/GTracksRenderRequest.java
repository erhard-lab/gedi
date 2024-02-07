package gedi.gui.gtracks.rendering;

public class GTracksRenderRequest {
	
	private double height;

	public GTracksRenderRequest(double height) {
		super();
		this.height = height;
	}

	public double getHeight() {
		return height;
	}

	@Override
	public String toString() {
		return "GTracksRenderRequest [height=" + height + "]";
	}
	
	

}
