package gedi.riboseq.visu;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.util.HashMap;
import java.util.function.Function;

import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.gui.genovis.style.StyleObject;
import gedi.riboseq.inference.orf.Orf;
import gedi.util.ArrayUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.gui.PixelBasepairMapper;

public class OrfRenderer extends PeptideRenderer<Orf> {

	
	
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
			Orf data, double xOffset, double y, double h, boolean boxes, boolean lines) {

		frameColors = brighterframeColors;
		if (data.hasStart() && data.hasStop())
			stringer = o->"";
		else
			stringer = o->o.getData().getOrfType().toString();
		
		Function<ReferenceGenomicRegion<Orf>, Paint> oldBorder = border;
		if (!data.passesAllFilters()) 
			border = (rgr)->Color.red;
		
		GenomicRegion re = super.renderBox(g2, locationMapper, reference, strand, region, data, xOffset, y, h, boxes, lines);
		
		stringer = o->o.getData().getOrfType().toString();
		frameColors = save;
		if (data.hasStart() && data.hasStop()) {
			GenomicRegion codingRegion = data.getStartToStop(reference.toStrand(strand), region);
			super.renderBox(g2, locationMapper, reference, strand, codingRegion, data, xOffset, y, h, boxes, lines);
		}
		
		border = oldBorder;
		
		return re;
	}

	
}
