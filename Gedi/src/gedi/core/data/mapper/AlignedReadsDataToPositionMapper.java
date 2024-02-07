package gedi.core.data.mapper;

import java.util.Iterator;
import java.util.Map.Entry;
import java.util.function.BiPredicate;

import javax.script.ScriptException;

import gedi.core.data.numeric.GenomicNumericProvider.SpecialAggregators;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.nashorn.JSBiPredicate;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=PixelBlockToValuesMap.class)
public class AlignedReadsDataToPositionMapper implements GenomicRegionDataMapper<IntervalTree<GenomicRegion,AlignedReadsData>, PixelBlockToValuesMap>{

	private GenomicRegionPosition referencePosition;
	private int offset;
	private SpecialAggregators aggregator = SpecialAggregators.Sum; 
	
	private ReadCountMode readCountMode = ReadCountMode.Weight;
	
	public AlignedReadsDataToPositionMapper() {
		this(GenomicRegionPosition.Start,0);
	}


	public AlignedReadsDataToPositionMapper(GenomicRegionPosition referencePosition) {
		this(referencePosition,0);
	}

	public AlignedReadsDataToPositionMapper(GenomicRegionPosition referencePosition,
			int offset) {
		this.referencePosition = referencePosition;
		this.offset = offset;
	}
	
	public void setReferencePosition(GenomicRegionPosition referencePosition) {
		this.referencePosition = referencePosition;
	}
	
	public void setReadCountMode(ReadCountMode readCountMode) {
		this.readCountMode = readCountMode;
	}

	public GenomicRegionPosition getReferencePosition() {
		return referencePosition;
	}

	public int getOffset() {
		return offset;
	}

	public void setAggregator(SpecialAggregators aggregator) {
		this.aggregator = aggregator;
	}
	
	private BiPredicate<AlignedReadsData,Integer> filter = null;
	public void setFilter(String js) throws ScriptException {
		this.filter = new JSBiPredicate<>(false, "function(data,d) "+js);
	}
	public void setSlamFilterSense(int num) {
		setHasMismatchFilter('T', 'C',num);
	}
	public void setSlamFilterAntisense(int num) {
		setHasMismatchFilter('A', 'G',num);
	}
	public void setHasMismatchFilter(char genomic, char read) {
		setHasMismatchFilter(genomic, read, 1);
	}
	public void setHasMismatchFilter(char genomic, char read, int num) {
		char igenomic = SequenceUtils.getDnaComplement(genomic);
		char iread = SequenceUtils.getDnaComplement(read);
		
		this.filter = (data,d)->{
			int c = 0;
			for (int v=0; v<data.getVariationCount(d); v++)
				if (data.isMismatch(d, v) && !data.isVariationFromSecondRead(d, v) && data.getMismatchGenomic(d, v).charAt(0)==genomic && data.getMismatchRead(d, v).charAt(0)==read)
					c++;
				else if (data.isMismatch(d, v) && data.isVariationFromSecondRead(d, v) && data.getMismatchGenomic(d, v).charAt(0)==igenomic && data.getMismatchRead(d, v).charAt(0)==iread)
					c++;
			return c>=num;
		};
	}

	

	@Override
	public PixelBlockToValuesMap map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			IntervalTree<GenomicRegion, AlignedReadsData> data) {
		if (data.isEmpty()) return new PixelBlockToValuesMap();
		
		int numCond = data.firstEntry().getValue().getNumConditions();

		PixelBlockToValuesMap re = new PixelBlockToValuesMap(pixelMapping, numCond, Double.NaN);
		MutableReferenceGenomicRegion rgr = new MutableReferenceGenomicRegion().setReference(data.getReference());
		int[] count = new int[re.size()];
		
		if (filter==null) {
			Iterator<Entry<GenomicRegion, AlignedReadsData>> it = data.entrySet().iterator();
			while (it.hasNext()) {
				Entry<GenomicRegion, AlignedReadsData> e = it.next();
				
				int block = pixelMapping.getBlockForBp(reference, referencePosition.position(rgr.setRegion(e.getKey()), offset));
				if (block==-1 || block>=re.size()) continue;
				
				NumericArray vals = re.getValues(block);
				AlignedReadsData ard = e.getValue();
				for (int c=0; c<numCond; c++) {
					double p = vals.getDouble(c);
					double rc = ard.getTotalCountForCondition(c, readCountMode);
					if (Double.isNaN(p))
						vals.set(c, rc);
					else
						vals.set(c, aggregator.increment(p, rc));
				}
				count[block]++;

			}
		} 
		else {
			Iterator<Entry<GenomicRegion, AlignedReadsData>> it = data.entrySet().iterator();
			while (it.hasNext()) {
				Entry<GenomicRegion, AlignedReadsData> e = it.next();
				
				int block = pixelMapping.getBlockForBp(reference, referencePosition.position(rgr.setRegion(e.getKey()), offset));
				if (block==-1 || block>=re.size()) continue;
				
				NumericArray vals = re.getValues(block);
				AlignedReadsData ard = e.getValue();
				
				for (int d=0; d<ard.getDistinctSequences(); d++) {
					if (filter.test(ard,d)) {
						for (int c=0; c<numCond; c++) {
							double p = vals.getDouble(c);
							double rc = ard.getCount(d,c,readCountMode);
							if (Double.isNaN(p))
								vals.set(c, rc);
							else
								vals.set(c, aggregator.increment(p, rc));
						}
					}
				}
					
				
				count[block]++;

			}
		}
		
		for (int i=0; i<pixelMapping.size(); i++) {
			NumericArray vals = re.getValues(i);
			for (int j=0; j<numCond; j++)
				vals.set(j, aggregator.getIncrementalResult(vals.getDouble(j), count[i]));
		}



		return re;
	}


	//
	//	@Override
	//	public GenomicNumericProvider map(ReferenceSequence reference,
	//			GenomicRegion region,PixelLocationMapping pixelMapping,
	//			IntervalTree<GenomicRegion, AlignedReadsData> data) {
	//		if (data.isEmpty()) return new EmptyGenomicNumericProvider();
	//		
	//		int numCond=data.firstEntry().getValue().getNumConditions();
	//		double dens = data.size()/(double)region.getTotalLength();
	//		
	//		if (dens<MAX_DENSE) {
	//			TreeMap<Integer,NumericArray> re = new TreeMap<Integer,NumericArray>();
	//			data.getIntervalsIntersecting(region, e->{
	//				int pos = referencePosition.position(reference, e.getKey(),offset);
	//				if (region.contains(pos)) {
	//					NumericArray pr = re.get(pos);
	//					if (pr==null) re.put(pos, pr = NumericArray.createMemory(numCond,NumericArrayType.Integer));
	//					e.getValue().addTotalCount(pr);
	//				}
	//			});
	//			return new SparseGenomicNumericProvider(reference, region, re);
	//		}
	//		
	//		NumericArray[] re = new NumericArray[numCond];
	//		for (int i=0; i<re.length; i++)
	//			re[i] = NumericArray.createMemory(region.getTotalLength(),NumericArrayType.Integer);
	//		
	//		data.getIntervalsIntersecting(region, e->{
	//			int pos = referencePosition.position(reference, e.getKey(),offset);
	//			if (region.contains(pos)) {
	//				int ind = region.induce(pos);
	//				for (int i=0; i<re.length; i++)
	//					re[i].add(ind, e.getValue().getTotalCount(i));
	//			}
	//		});
	//		
	//		return new DenseGenomicNumericProvider(reference,region,re);
	//	}


}
