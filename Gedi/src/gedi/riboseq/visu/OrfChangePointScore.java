package gedi.riboseq.visu;


import gedi.core.data.mapper.GenomicRegionDataMapping;
import gedi.core.reference.Strand;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.riboseq.inference.orf.Orf;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=PixelBlockToValuesMap.class)
public class OrfChangePointScore extends OrfToScore {

	public OrfChangePointScore(Strand strand) {
		super(strand);
 	}

	@Override
	protected double getScore(Orf orf, int p) {
		return orf.getChangePointScores()[p];
	}


}
