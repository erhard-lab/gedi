package gedi.core.data.mapper;

import java.util.Iterator;
import java.util.Map.Entry;

import gedi.core.data.numeric.GenomicNumericProvider.SpecialAggregators;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.sequence.Alphabet;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=PixelBlockToValuesMap.class)
public class AlignedReadsDataToSoftclipMapper implements GenomicRegionDataMapper<IntervalTree<GenomicRegion,AlignedReadsData>, PixelBlockToValuesMap>{

	private ReadCountMode mode = ReadCountMode.Weight;
	
	private SpecialAggregators aggregator = SpecialAggregators.Sum; 


	public void setAggregator(SpecialAggregators aggregator) {
		this.aggregator = aggregator;
	}
	
	public void setReadCountMode(ReadCountMode mode) {
		this.mode = mode;
	}

	private boolean p3=true;
	public void setFivePrime() {
		p3=false;
	}
	
	
	@Override
	public PixelBlockToValuesMap map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			IntervalTree<GenomicRegion, AlignedReadsData> data) {
		if (data.isEmpty()) return new PixelBlockToValuesMap();
		

		PixelBlockToValuesMap re = new PixelBlockToValuesMap(pixelMapping, 4, 0);
		MutableReferenceGenomicRegion rgr = new MutableReferenceGenomicRegion().setReference(data.getReference());
		int[] count = new int[re.size()];
		
		int[] ind = Alphabet.getDna().createIndex();
		
		Iterator<Entry<GenomicRegion, AlignedReadsData>> it = data.entrySet().iterator();
		while (it.hasNext()) {
			Entry<GenomicRegion, AlignedReadsData> e = it.next();
			
			AlignedReadsData ard = e.getValue();
			rgr.setRegion(e.getKey());
					
			int pos = rgr.getReference().isMinus()==p3?rgr.getRegion().getStart():rgr.getRegion().getStop();
			double[] avals = new double[4];
			
			int block = pixelMapping.getBlockForBp(reference, pos);
			if (block<0 || block>=pixelMapping.size()) continue;

			NumericArray vals = re.getValues(block);
			
			
			for (int d=0; d<ard.getDistinctSequences(); d++) {
				
				double rc = ard.getTotalCountForDistinct(d, mode);
				for (int v=0; v<ard.getVariationCount(d); v++) {
					
					if (ard.isSoftclip(d, v) && checkSide(ard,d,v)) {
						CharSequence seq = ard.getSoftclip(d, v);
						
						for (int sp=0; sp<seq.length(); sp++) {
							int t = ind[seq.charAt(sp)];
							if (t>=0 && t<4) {
								avals[t]+=rc;
							}
						}
					}
				}
			}
			
			for (int t=0; t<4; t++) {
				double p = vals.getDouble(t);
				if (Double.isNaN(p))
					vals.set(t, avals[t]);
				else
					vals.set(t, aggregator.increment(p, avals[t]));
			}
			count[block]++;
						

		}
		
		for (int i=0; i<pixelMapping.size(); i++) {
			NumericArray vals = re.getValues(i);
			for (int j=0; j<vals.length(); j++) {
				double res = aggregator.getIncrementalResult(vals.getDouble(j), count[i]);
				vals.set(j, res);
			}
		}



		return re;
	}

	private boolean checkSide(AlignedReadsData ard, int d, int v) {
		if (p3) return !ard.isSoftclip5p(d, v);
		else return ard.isSoftclip5p(d, v);
	}



}
