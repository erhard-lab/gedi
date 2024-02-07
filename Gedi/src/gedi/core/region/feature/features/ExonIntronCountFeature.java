package gedi.core.region.feature.features;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.ReferenceSequenceConversion;
import gedi.core.reference.Strand;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Set;
import java.util.function.BiPredicate;
import java.util.function.UnaryOperator;


@GenomicRegionFeatureDescription(toType=Object.class)
public class ExonIntronCountFeature<T> extends AbstractFeature<Object> {

	private Genomic genomic;
	private ReferenceSequenceConversion referenceSequenceConversion = ReferenceSequenceConversion.none;
	
	
	@Override
	public GenomicRegionFeature<Object> copy() {
		ExonIntronCountFeature<T> re = new ExonIntronCountFeature<T>(genomic);
		re.copyProperties(this);
		re.referenceSequenceConversion = referenceSequenceConversion;
		return re;
	}
	
	public ExonIntronCountFeature(Genomic genomic) {
		this.genomic = genomic;
		minValues = 1;
		maxValues = 1;
		minInputs = maxInputs = 0;
	}
	
	
	public void setReferenceSequenceConversion(
			ReferenceSequenceConversion referenceSequenceConversion) {
		this.referenceSequenceConversion = referenceSequenceConversion;
	}
	

	@Override
	protected void accept_internal(Set<Object> values) {
		ReferenceSequence reference = referenceSequenceConversion.apply(this.referenceRegion.getReference());
		
		boolean exon = false;
		boolean gene = false;
		
		AlignedReadsData ard = (AlignedReadsData) referenceRegion.getData();
		
		
		gene = genomic.getGenes().ei(reference, referenceRegion.getRegion())
			.filter(rgr->rgr.getRegion().contains(referenceRegion.getRegion()))
			.count()>0;
		
		
		exon = genomic.getTranscripts().ei(reference, referenceRegion.getRegion())
			.filter(rgr->ard.isConsistentlyContained(referenceRegion, rgr, 0))
			.count()>0;
				
			
		
		if (reference.getStrand()==Strand.Independent) {
			gene |= genomic.getGenes().ei(reference.toPlusStrand(), referenceRegion.getRegion())
					.chain(genomic.getGenes().ei(reference.toMinusStrand(), referenceRegion.getRegion()))
					.filter(rgr->rgr.getRegion().contains(referenceRegion.getRegion()))
					.count()>0;
					
			exon |= genomic.getTranscripts().ei(reference.toPlusStrand(), referenceRegion.getRegion())
				.chain(genomic.getTranscripts().ei(reference.toMinusStrand(), referenceRegion.getRegion()))
				.filter(rgr->ard.isConsistentlyContained(referenceRegion, rgr, 0))
				.count()>0;
		}
		
		
		if (exon) values.add("Exonic");
		else if (gene) values.add("Intronic");
		else values.add("Intergenic");
		
	}
	
}
