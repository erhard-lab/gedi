package gedi.gui.genovis.tracks.boxrenderer;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;

import gedi.core.data.annotation.NarrowPeakAnnotation;
import gedi.core.data.annotation.ScoreProvider;
import gedi.core.data.mapper.StorageSource;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.GenomicRegion;
import gedi.util.MathUtils;
import gedi.util.PaintUtils;
import gedi.util.functions.EI;
import gedi.util.gui.PixelBasepairMapper;
import gedi.util.gui.ValueToColorMapper;

public class NarrowPeakRenderer extends BoxRenderer<NarrowPeakAnnotation> {

	private ValueToColorMapper mapper = new ValueToColorMapper(Color.WHITE,Color.BLACK);
	
	private Color summitColor=Color.RED;
	
	public NarrowPeakRenderer() {
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
	
	public void summit(Color color) {
		this.summitColor=color;
	}


	@Override
	public GenomicRegion renderBox(Graphics2D g2, PixelBasepairMapper locationMapper,
			ReferenceSequence reference, Strand strand, GenomicRegion region,
			NarrowPeakAnnotation d, double xOffset, double y, double h,boolean boxes, boolean lines) {
		GenomicRegion re = super.renderBox(g2, locationMapper, reference, strand, region, d, xOffset, y, h,boxes,lines);
		
		if (d.getSummit()>=0 && d.getSummit()<region.getTotalLength()) {
			Rectangle2D tile = getTile(reference, region.map(d.getSummit()), region.map(d.getSummit())+1, locationMapper, xOffset, y, h);
			g2.setPaint(summitColor);
			g2.fill(tile);
			g2.draw(tile);
		}
		
		return re;
	}

}
