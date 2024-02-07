package gedi.gui.genovis.tracks.boxrenderer;

import gedi.util.PaintUtils;

import java.awt.Color;

public class AnnotationRenderer<T> extends BoxRenderer<T> {

	
	
	public AnnotationRenderer() {
		setHeight(20);
		setFont("Arial", 14, true, false);
		setForeground(t->Color.WHITE);
		setBackground(t->PaintUtils.parseColor("#004d6dff"));
		stringer = x->String.valueOf(x.getData());
	}
	


}
