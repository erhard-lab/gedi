package gedi.core.region.feature.features;

import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;

import java.util.Set;


@GenomicRegionFeatureDescription(fromType=Void.class,toType=String.class)
public class GenomeFeature extends AbstractFeature<String> {

	private Genomic g;
	
	public GenomeFeature(Genomic g) {
		minInputs = maxInputs = 0;
		this.g = g;
	}
	
	@Override
	public GenomicRegionFeature<String> copy() {
		GenomeFeature re = new GenomeFeature(g);
		re.copyProperties(this);
		return re;
	}


	@Override
	protected void accept_internal(Set<String> values) {
		
		Genomic or = g.getOrigin(referenceRegion.getReference());
		if (or!=null)
			values.add(or.getId());
		
	}
	

}
