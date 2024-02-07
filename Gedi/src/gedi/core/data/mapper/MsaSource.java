package gedi.core.data.mapper;

import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.sequence.CompositeSequenceProvider;
import gedi.core.sequence.FastaIndexSequenceProvider;
import gedi.core.sequence.SequenceProvider;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.algorithm.string.alignment.multiple.MsaBlock;
import gedi.util.algorithm.string.alignment.multiple.MultipleSequenceAlignment;
import gedi.util.dynamic.DynamicObject;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.function.Consumer;

@GenomicRegionDataMapping(fromType=Void.class,toType=CharSequence[].class)
public class MsaSource implements GenomicRegionDataSource<CharSequence[]>{
	
	private MultipleSequenceAlignment msa;
	
	public MsaSource(GenomicRegionStorage<MsaBlock> storage) {
		this.msa = new MultipleSequenceAlignment(storage);
	}
	
	@Override
	public CharSequence[] get(ReferenceSequence reference, GenomicRegion region,PixelLocationMapping pixelMapping) {
		return msa.getMsa(new ImmutableReferenceGenomicRegion<>(reference, region));
	}

	private String id = null;
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}


	@Override
	public DynamicObject getMeta() {
		LinkedHashMap<String, Object> meta = new LinkedHashMap<>();
		meta.put("species", msa.getSpecies());
		return DynamicObject.from(meta);
	}
	
	@Override
	public <T> void applyForAll(Class<T> cls, Consumer<T> consumer) {
	}
	
}
