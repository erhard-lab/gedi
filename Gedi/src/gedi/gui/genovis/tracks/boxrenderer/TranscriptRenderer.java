package gedi.gui.genovis.tracks.boxrenderer;

import java.awt.Color;
import java.awt.Graphics2D;

import gedi.core.data.annotation.Transcript;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.GenomicRegion;
import gedi.util.PaintUtils;
import gedi.util.gui.PixelBasepairMapper;

public class TranscriptRenderer extends BoxRenderer<Transcript> {

	
	private BoxRenderer<Transcript> coding = new BoxRenderer<Transcript>();
	
	public TranscriptRenderer() {
		coding.setBackground(t->PaintUtils.parseColor("#5b95ad"));
		coding.setForeground(t->Color.WHITE);
		coding.setHeight(20);
		coding.setFont("Arial", 14, true, false);
		coding.stringer = t->t.getData().getTranscriptId();
		setHeight(20);
		setFont("Arial", 14, true, false);
		setForeground(t->Color.WHITE);
		setBackground(t->PaintUtils.parseColor("#98c1d2"));
		stringer = d->"";
	}
	

	@Override
	public GenomicRegion  renderBox(Graphics2D g2,
			PixelBasepairMapper locationMapper, ReferenceSequence reference, Strand strand, GenomicRegion region,
			Transcript data, double xOffset, double y, double h,boolean boxes, boolean lines) {
		stringer = data.isCoding()?null:coding.stringer;
			
		GenomicRegion re = super.renderBox(g2, locationMapper, reference, strand, region, data, xOffset, y, h,boxes,lines);
		if (data.isCoding()) {
			GenomicRegion codingRegion = data.getCds(reference, region);
			coding.renderBox(g2, locationMapper, reference, strand, codingRegion, data, xOffset, y, h,boxes,lines);
		}
		
		return re;
	}


}
