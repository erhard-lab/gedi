package gedi.gui.genovis.tracks.boxrenderer;

import gedi.core.data.annotation.ScoreProvider;
import gedi.core.data.mapper.StorageSource;
import gedi.util.MathUtils;
import gedi.util.PaintUtils;
import gedi.util.functions.EI;
import gedi.util.gui.ColorPalettes;
import gedi.util.gui.ValueToColorMapper;

import java.awt.Color;

public class ScoreRenderer extends BoxRenderer<ScoreProvider> {

	private ValueToColorMapper mapper = new ValueToColorMapper(Color.WHITE,Color.BLACK);
	
	public ScoreRenderer() {
		setBackground(t->mapper.apply(t.getScore()));
		setForeground(t->PaintUtils.isDarkColor(mapper.apply(t.getScore()))?Color.WHITE:Color.BLACK);
		stringer = s->s.getData().toString();
	}
	
	public void linear(double min, double max) {
		mapper = new ValueToColorMapper(MathUtils.linearRange(min,max), mapper.getColors());
	}
	
	public void linear(StorageSource<? extends ScoreProvider> src) {
		double[] minmax = EI.wrap(src.getStorages()).demultiplex(s->s.ei()).mapToDouble(r->r.getData().getScore()).saveMinMax(0,1000);
		linear(minmax[0],minmax[1]);
	}
	
	
	public void colors(Color... colors) {
		mapper = new ValueToColorMapper(mapper.getRange(), colors);
	}
	
	public void colors(double min, double max, Color... colors) {
		mapper = new ValueToColorMapper(MathUtils.linearRange(min, max), colors);
	}
	
	public void colors(ColorPalettes palette, int numColors) {
		colors(palette.getPalette(numColors));
	}
	

}
