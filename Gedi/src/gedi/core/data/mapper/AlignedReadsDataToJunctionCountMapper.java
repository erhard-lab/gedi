package gedi.core.data.mapper;


import java.util.function.BiPredicate;

import javax.script.ScriptException;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.MissingInformationIntronInformation;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.nashorn.JSBiPredicate;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class AlignedReadsDataToJunctionCountMapper implements GenomicRegionDataMapper<IntervalTree<GenomicRegion,AlignedReadsData>, IntervalTree<GenomicRegion,NumericArray>>{

	private Strand strand = Strand.Independent;
	private ReadCountMode readCountMode = ReadCountMode.Weight;
	
	public void setReadCountMode(ReadCountMode readCountMode) {
		this.readCountMode = readCountMode;
	}
	
	public void setStrand(Strand strand) {
		this.strand = strand;
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
		setHasMismatchFilter(genomic,read,1);
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
	public IntervalTree<GenomicRegion,NumericArray> map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			IntervalTree<GenomicRegion, AlignedReadsData> data) {
		if (data.isEmpty()) return new IntervalTree<GenomicRegion,NumericArray>(data.getReference());
		
		int numCond = data.firstEntry().getValue().getNumConditions();
		
		IntervalTree<GenomicRegion,NumericArray> re = new IntervalTree<GenomicRegion,NumericArray>(data.getReference());

		if (filter==null)
			data.entrySet().iterator().forEachRemaining(e->{
				if (e.getKey().getNumParts()>1) {
					int l=!strand.isMinus()?0:e.getKey().getTotalLength();
					for (int i=0; i<e.getKey().getNumParts()-1; i++) {
						if (!strand.isMinus())
							l+=e.getKey().getLength(i);
						else
							l-=e.getKey().getLength(i);
						if (e.getKey() instanceof MissingInformationIntronInformation && ((MissingInformationIntronInformation)e.getKey()).isMissingInformationIntron(i)){}
						else {
							AlignedReadsData val = e.getValue();
							if (!val.isFalseIntron(l,0)) {
								ArrayGenomicRegion reg = new ArrayGenomicRegion(e.getKey().getEnd(i),e.getKey().getStart(i+1));
								NumericArray a = re.computeIfAbsent(reg, k->NumericArray.createMemory(numCond, NumericArrayType.Double));
								val.addTotalCountsForConditions(a, readCountMode);
							}
						}
					}
				}
			});
		else
			data.entrySet().iterator().forEachRemaining(e->{
				if (e.getKey().getNumParts()>1) {
					AlignedReadsData ard = e.getValue();
					for (int d=0; d<ard.getDistinctSequences(); d++) {
						if (filter.test(ard,d)) {
							int l=!strand.isMinus()?0:e.getKey().getTotalLength();
							for (int i=0; i<e.getKey().getNumParts()-1; i++) {
								if (!strand.isMinus())
									l+=e.getKey().getLength(i);
								else
									l-=e.getKey().getLength(i);
								if (e.getKey() instanceof MissingInformationIntronInformation && ((MissingInformationIntronInformation)e.getKey()).isMissingInformationIntron(i)){}
								else {
									
									if (!ard.isFalseIntron(l,0)) {
										ArrayGenomicRegion reg = new ArrayGenomicRegion(e.getKey().getEnd(i),e.getKey().getStart(i+1));
										NumericArray a = re.computeIfAbsent(reg, k->NumericArray.createMemory(numCond, NumericArrayType.Double));
										ard.addCountsForDistinct(d, a, readCountMode);
									}
								}
							}
						}
					}
				}
			});
		return re;
	}



}
