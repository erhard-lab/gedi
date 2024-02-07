package gedi.riboseq.visu;


import gedi.core.data.mapper.GenomicRegionDataMapping;
import gedi.core.reference.Strand;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.riboseq.inference.orf.Orf;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=PixelBlockToValuesMap.class)
public class PriceOrfChangePoint extends PriceOrfToScore {

	int cond;
	public PriceOrfChangePoint(Strand strand, double cond) {
		super(strand);
		this.cond = (int)cond;
	}

	@Override
	protected double getScore(PriceOrf orf, int p) {
		double psum = 0;
		double asum = 0;
		for (int x=p+10; x>=p-10; x--) {
			if (x>=0 && x<orf.getNumPositions())
				asum = asum+orf.getProfile(cond,x)+0.1;
			if (x==p) psum = asum;
		}
		return psum/asum;
	}

	@Override
	protected int getRows(PriceOrf orf) {
		return 1;
	}

}
