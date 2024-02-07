package gedi.util.genomic.metagene;

import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;

public class ExonIntronMetageneDataProvider<D> implements MetageneDataProvider {

	private MemoryIntervalTreeStorage<?> transcripts;
	private Strandness strandness = Strandness.Sense;
	

	public ExonIntronMetageneDataProvider(Genomic genomic) {
		this.transcripts = genomic.getTranscripts();
	}
	
	public ExonIntronMetageneDataProvider(GenomicRegionStorage<?> transcripts) {
		this.transcripts = transcripts.toMemory();
	}
	
	public ExonIntronMetageneDataProvider<D> setStrandness(Strandness strandness) {
		this.strandness = strandness;
		return this;
	}
	

	@Override
	public AutoSparseDenseDoubleArrayCollector getData(ImmutableReferenceGenomicRegion<?> rgr) {
		AutoSparseDenseDoubleArrayCollector re = new AutoSparseDenseDoubleArrayCollector(rgr.getRegion().getTotalLength()/10,rgr.getRegion().getTotalLength());
		
		if (!strandness.equals(Strandness.Antisense))
			for (ImmutableReferenceGenomicRegion<?> tr : transcripts.ei(rgr).loop()) {
				for (int p=0; p<tr.getRegion().getNumParts(); p++) {
					for (int i=tr.getRegion().getStart(p); i<tr.getRegion().getEnd(p); i++)
						if (rgr.getRegion().contains(i))
							re.add(rgr.induce(i),1);
				}
			}
		else if (!strandness.equals(Strandness.Sense))
			for (ImmutableReferenceGenomicRegion<?> tr : transcripts.ei(rgr.toOppositeStrand()).map(t->t.toOppositeStrand()).loop()) {
				for (int p=0; p<tr.getRegion().getNumParts(); p++) {
					for (int i=tr.getRegion().getStart(p); i<tr.getRegion().getEnd(p); i++)
						if (rgr.getRegion().contains(i))
							re.add(rgr.induce(i),1);
				}
			}

		return re;
	}

	
	
}
