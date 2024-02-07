package gedi.gui.genovis.tracks.boxrenderer;

import java.awt.Color;

import gedi.gui.genovis.VisualizationTrack;

public class TrackLabelRenderer<T> extends BoxRenderer<T> {

	
	
	public TrackLabelRenderer(VisualizationTrack<?, ?> track) {
		setForeground(t->Color.WHITE);
		setBackground(t->Color.BLACK);
		stringer = x->track.getId();
	}
	
	public TrackLabelRenderer() {
		setForeground(t->Color.WHITE);
		setBackground(t->Color.BLACK);
	}


}
