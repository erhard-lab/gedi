package gedi.core.genomic;

import gedi.core.data.annotation.Transcript;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.datastructure.tree.Trie;
import java.util.function.Function;

public class NameIndex {

	private Genomic genomic;
	private Trie<ReferenceGenomicRegion<?>> trie;
	
	public NameIndex(String a) {} // dummy constructor for compatibility reasons
	
	
	public NameIndex(Genomic g) {
		this.genomic = g;
	}
	
	
	public Trie<ReferenceGenomicRegion<?>> getIndex() {
		if (trie==null)
			synchronized (this) {
				if (trie==null) {
					trie = new Trie<>();
					Function<String,String> fun = genomic.getGeneTableColumns().contains("symbol")?genomic.getGeneTable("symbol"):null;
					for (ImmutableReferenceGenomicRegion<String> g : genomic.getGenes().ei().loop()) {
						trie.put(g.getData(), g);
						if (fun!=null && fun.apply(g.getData())!=null)
							trie.put(fun.apply(g.getData()), g);
					}
					for (ImmutableReferenceGenomicRegion<Transcript> t : genomic.getTranscripts().ei().loop()) {
						trie.put(t.getData().getTranscriptId(), t);
					}
				}
			}
		return trie;
	}


	
}
