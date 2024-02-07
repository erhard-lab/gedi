package gedi.core.data.mapper;

import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

import java.util.function.Function;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class GenomicStorageFilter<D> extends StorageFilter<D >{

	private Genomic genomic;
	
	public GenomicStorageFilter(Genomic genomic) {
		this.genomic = genomic;
	}

	public void addTranscriptTable(String col) {
		addFunction(genomic.getTranscriptTable(col));
	}
	
	public void addGeneTable(String col) {
		addFunction(genomic.getGeneTable(col));
	}
	
	public void addFunction(Function<String, String> map) {
		if (map!=null)
		addOperator(rgr->rgr.toMutable().transformData(d->{
			if (d==null) return (D)map.apply("");
			if (d instanceof Transcript) return (D)map.apply(((Transcript)d).getTranscriptId());
			return (D)map.apply(d.toString());
		}));		
	}

}
