package gedi.core.data.mapper;

import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.sequence.CompositeSequenceProvider;
import gedi.core.sequence.FastaIndexSequenceProvider;
import gedi.core.sequence.SequenceProvider;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;

import java.io.IOException;
import java.util.function.Consumer;

@GenomicRegionDataMapping(fromType=Void.class,toType=CharSequence.class)
public class SequenceSource implements GenomicRegionDataSource<CharSequence>{
	
	private CompositeSequenceProvider prov = new CompositeSequenceProvider();
	
	public SequenceSource() {
	}

	public void addFastaIndex(String file) throws IOException {
		prov.add(new FastaIndexSequenceProvider(file));
	}
	
	public void addGenomic(Genomic g) throws IOException {
		prov.add(g);
	}
	
	public void addName(String name) throws IOException {
		prov.add(Genomic.get(name));
	}
	
	@Override
	public CharSequence get(ReferenceSequence reference, GenomicRegion region,PixelLocationMapping pixelMapping) {
		return prov.getSequence(reference, region);
	}

	private String id = null;
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}

	
	@Override
	public <T> void applyForAll(Class<T> cls, Consumer<T> consumer) {
		for (SequenceProvider p : prov.getProviders())
			if (cls.isInstance(p))
				consumer.accept(cls.cast(p));
	}
	
}
