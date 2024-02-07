package gedi.riboseq.analysis;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.features.AbstractFeature;
import gedi.riboseq.inference.orf.Orf;
import gedi.riboseq.inference.orf.PriceOrf;

import java.util.Set;


@GenomicRegionFeatureDescription(toType=String.class)
public class OrfTypeFeature extends AbstractFeature<String> {

	public OrfTypeFeature() {
		minValues = maxValues = 1;
		minInputs = maxInputs = 0;
	}

	@Override
	protected void accept_internal(Set<String> values) {
		PriceOrf orf = (PriceOrf) referenceRegion.getData();
		values.add(orf.getType().name());
	}

	@Override
	public GenomicRegionFeature<String> copy() {
		OrfTypeFeature re = new OrfTypeFeature();
		re.copyProperties(this);
		return re;
	}
	
}
