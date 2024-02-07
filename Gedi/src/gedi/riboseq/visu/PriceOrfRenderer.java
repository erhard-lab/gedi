package gedi.riboseq.visu;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.util.HashMap;
import java.util.function.Function;

import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.gui.genovis.style.StyleObject;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.util.ArrayUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.gui.PixelBasepairMapper;

public class PriceOrfRenderer extends PeptideRenderer<PriceOrf> {

	
	
	private Color[] white = {Color.WHITE,Color.WHITE,Color.WHITE};
	private Color[] black = {Color.black,Color.black,Color.black};
	private Color[] fnull = {null,null,null};
	private Color[] gray = {Color.gray,Color.gray,Color.gray};
	private Color[] brighterframeColors = {Color.WHITE,Color.WHITE,Color.WHITE};
	private Color[] save = {Color.WHITE,Color.WHITE,Color.WHITE};


	
	public void setStyles(DynamicObject styles) {
		if (styles.isArray()){
			HashMap<String, Color> colorMap = ArrayUtils.createMapping(styles.applyTo(new StyleObject[styles.length()]),s->s.getName(),s->s.getColor());
			for (int f=0; f<3;f++)
				brighterframeColors[f] = colorMap.getOrDefault("Frame"+f,Color.white).brighter();
		}
		super.setStyles(styles);
		save = frameColors;
	}
	
	
	@Override
	public GenomicRegion renderBox(Graphics2D g2,
			PixelBasepairMapper locationMapper, ReferenceSequence reference, Strand strand, GenomicRegion region,
			PriceOrf data, double xOffset, double y, double h, boolean boxes, boolean lines) {

		stringer = o->"";
		frameColors = save;
		GenomicRegion codingRegion = data.getStartStop(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand),region),0,true).getRegion();
		super.renderBox(g2, locationMapper, reference, strand, codingRegion, data, xOffset, y, h, boxes, lines);

		
//		frameColors = white;
//		super.renderBox(g2, locationMapper, reference, strand, region, data, xOffset, y, h);
		
//		TriFunction<ReferenceSequence, GenomicRegion, PriceOrf, Paint> oldBorder = border;
		
//		frameColors = brighterframeColors;
		for (int a=0; a<data.getNumAlternativeStartCodons(); a++) {
			
			codingRegion = data.getStartStop(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand),region),a,true).getRegion();
			GenomicRegion atg = new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand),codingRegion).map(new ArrayGenomicRegion(0,3));
			
			if (data.isPredictedStartIndex(a)) 
				frameColors = black;
			else
				frameColors = gray;
			super.renderBox(g2, locationMapper, reference, strand, atg, data, xOffset, y, h, boxes, lines);
			
//			}
		}
		
		Function<ReferenceGenomicRegion<PriceOrf>, Paint> borderSave = this.border;
		Function<ReferenceGenomicRegion<PriceOrf>, Paint> backgroundSave = this.background;
		
		this.border = null;
		frameColors = fnull;
		stringer = o->o.getData().getTranscript();
		codingRegion = data.getStartStop(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand),region),0,true).getRegion();
		super.renderBox(g2, locationMapper, reference, strand, codingRegion, data, xOffset, y, h, boxes, lines);
		
		this.border = borderSave;
		
		frameColors = save;
		
//		setStringer(o->o.getTranscript());
//		GenomicRegion codingRegion = data.getStartStop(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand),region),true).getRegion();
//		super.renderBox(g2, locationMapper, reference, strand, codingRegion, data, xOffset, y, h);
		
//		border = oldBorder;
		
		return region;
	}

	
}
