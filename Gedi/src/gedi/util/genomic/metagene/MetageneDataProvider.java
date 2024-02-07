package gedi.util.genomic.metagene;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.datastructure.array.IndexDoubleProcessor;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;

public interface MetageneDataProvider {

	/**
	 * Size of the returned array usually is 1, only for the fast coverage algorithm, it is 2!
	 * @param ref
	 * @return
	 */
	IndexDoubleProcessor getData(ImmutableReferenceGenomicRegion<?> ref);
	
	public static MetageneDataProvider getTestV() {
		return (ref)->{
			AutoSparseDenseDoubleArrayCollector re = new AutoSparseDenseDoubleArrayCollector(ref.getRegion().getTotalLength(),ref.getRegion().getTotalLength());
			for (int i=0; i<ref.getRegion().getTotalLength(); i++) 
				re.add(i, Math.abs(1.0-i/(ref.getRegion().getTotalLength()-1.0)*2));
			return re;
		};
	}
	
}
