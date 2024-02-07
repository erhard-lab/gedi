package gedi.gui.genovis.tracks.boxrenderer;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataDecorator;
import gedi.core.data.reads.HasSubreads;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.gui.PixelBasepairMapper;
import gedi.util.gui.ValueToColorMapper;

public class AlignedReadsDataBoxRenderer extends BoxRenderer<AlignedReadsData> {

	private double factor = 1;
	private int[] conditions;
	
	private int subreads = 3;
	private Color[] subreadColors = createSubreadColors();
	
	private ReadCountMode readCountMode = ReadCountMode.Weight;
	
	public AlignedReadsDataBoxRenderer() {
		setBorder();
		setBackground(r->Color.BLACK);
	}
	
	private Color[] createSubreadColors() {
		ValueToColorMapper m = new ValueToColorMapper(Color.GRAY, Color.BLACK);
		Color[] re = new Color[subreads];
		for (int i=0; i<re.length; i++)
			re[i] = m.apply(1.0/(re.length-1)*i);
		return re;
	}

	public void setConditions(int[] conditions) {
		this.conditions = conditions;
	}
	
	public void setFactor(double factor) {
		this.factor = factor;
	}
	
	public void setSubreads(int subreads) {
		this.subreads = subreads;
		subreadColors = createSubreadColors();
	}
	
	public void setReadCountMode(ReadCountMode readCountMode) {
		this.readCountMode = readCountMode;
	}
	
	public ReadCountMode getReadCountMode() {
		return readCountMode;
	}

	private HasSubreads getSubreads(AlignedReadsData d) {
		if (d instanceof HasSubreads)
			return (HasSubreads) d;
		else if (d instanceof AlignedReadsDataDecorator)
			return getSubreads(((AlignedReadsDataDecorator) d).getParent());
		return null;
	}
	
	@Override
	public GenomicRegion renderBox(Graphics2D g2,
			PixelBasepairMapper locationMapper, ReferenceSequence reference, Strand strand, GenomicRegion region,
			AlignedReadsData data, double xOffset, double y, double h,boolean boxes, boolean lines) {
		

		GenomicRegion re;
		
		HasSubreads sr = getSubreads(data);
		
		if (sr!=null) {
			re = super.renderBox(g2, locationMapper, reference, strand, region, data, xOffset, y, h,false,true);
			
			ImmutableReferenceGenomicRegion rr = new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand), region);
			for (int s=0; s<sr.getNumSubreads(0); s++) {
				GenomicRegion srr = rr.map(new ArrayGenomicRegion(sr.getSubreadStart(0, s),sr.getSubreadEnd(0, s, region.getTotalLength())));
				
				
				setBackground(subreadColors[Math.min(subreads-1, sr.getSubreadId(0, s))]);
				super.renderBox(g2, locationMapper, reference, strand, srr, data, xOffset, y, h,true,false);
				
			}
		}
		else if (data.hasGeometry()) {
			
			GenomicRegion r1 = data.extractRead1(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand), region), 0).getRegion();
			GenomicRegion r2 = data.extractRead2(new ImmutableReferenceGenomicRegion<>(reference.toStrand(strand), region), 0).getRegion();
			
			setBackground(r->Color.GRAY);
			re = super.renderBox(g2, locationMapper, reference, strand, r1.union(r2), data, xOffset, y, h,true,true);
			setBackground(r->Color.BLACK);
			super.renderBox(g2, locationMapper, reference, strand, r1.intersect(r2), data, xOffset, y, h,true,true);
			
			
		} else
			re = super.renderBox(g2, locationMapper, reference, strand, region, data, xOffset, y, h,true,true);
		
		for (int d=0; d<data.getDistinctSequences(); d++) {
			double th = 0;
			
			
			double[] totals = data.getCountsForDistinct(d, readCountMode);//data.getTotalCountForDistinctSequence(d);
			
			if (conditions==null)
				th = ArrayUtils.sum(totals);
			else
				for (int c : conditions)
					th+=totals[c];
			
			th*=factor;
			
			if (th>0 && strand!=Strand.Independent) {
				for (int v=0; v<data.getVariationCount(d); v++) {
					if(data.isMismatch(d, v)) {
						double s = locationMapper.bpToPixel(reference,data.positionToGenomic(data.getMismatchPos(d, v),reference.toStrand(strand),region));
						double e = locationMapper.bpToPixel(reference,data.positionToGenomic(data.getMismatchPos(d, v),reference.toStrand(strand),region)+1);
						
						g2.setPaint(SequenceUtils.getNucleotideColorizer().apply(data.getMismatchRead(d, v).charAt(0)));
						Shape tile = new Rectangle2D.Double(xOffset+s, y, e-s, th-1);
						
						if (data.hasGeometry() && data.isPositionInOverlap(d, data.getMismatchPos(d, v))) {
							if (data.isVariationFromSecondRead(d, v))
								tile = getLowerTriangle((Rectangle2D)tile);
							else
								tile = getUpperTriangle((Rectangle2D)tile);
						}
						g2.draw(tile);
						g2.fill(tile);
					}
				}
				y+=th;
			}
		}
		
		return re;
	}

	private Shape getUpperTriangle(Rectangle2D tile) {
		GeneralPath re = new GeneralPath();
		re.moveTo(tile.getX(), tile.getY());
		re.lineTo(tile.getMaxX(), tile.getY());
		re.lineTo(tile.getX(), tile.getMaxY());
		re.closePath();
		return re;
	}

	private Shape getLowerTriangle(Rectangle2D tile) {
		GeneralPath re = new GeneralPath();
		re.moveTo(tile.getMaxX(), tile.getY());
		re.lineTo(tile.getMaxX(), tile.getMaxY());
		re.lineTo(tile.getX(), tile.getMaxY());
		re.closePath();
		return re;
	}

	@Override
	public double prefHeight(ReferenceSequence ref, GenomicRegion reg, AlignedReadsData d) {
		double re = 0;
		if (conditions==null)
			for (int i=0; i<d.getNumConditions(); i++)
				re+=d.getTotalCountForCondition(i, readCountMode);
		else 
			for (int i : conditions)
				re+=d.getTotalCountForCondition(i, readCountMode);
		return re*factor;
	}

}
