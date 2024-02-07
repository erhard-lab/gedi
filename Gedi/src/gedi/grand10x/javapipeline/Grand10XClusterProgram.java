package gedi.grand10x.javapipeline;

import java.util.ArrayList;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;

import executables.Grand10X;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.grand10x.TrimmedGenomicRegion;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.functions.EI;
import gedi.util.functions.IterateIntoSink;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand10XClusterProgram extends GediProgram {



	public Grand10XClusterProgram(Grand10xParameterSet params) {
		addInput(params.prefix);
		addInput(params.reads);
		addInput(params.genomic);
		addInput(params.nthreads);
		addInput(params.maxdist);

		addOutput(params.clustercit);
	}

	@Override
	public String execute(GediProgramContext context) throws Exception {
		String prefix = getParameter(0);
		GenomicRegionStorage<BarcodedAlignedReadsData> reads = getParameter(1);
		Genomic genomic = getParameter(2);
		int nthreads = getIntParameter(3);
		int maxdist = getIntParameter(4);
		
		genomic.getTranscriptTable("source");
		
		context.getLog().info("Identifying clusters...");
		
		TreeSet<ReferenceSequence> refs = new TreeSet<ReferenceSequence>(reads.getReferenceSequences());
		refs.removeIf(r->!genomic.getSequenceNames().contains(r.getName()));
		
		CenteredDiskIntervalTreeStorage<MutableInteger> out = new CenteredDiskIntervalTreeStorage<>(getOutputFile(0).getPath(),MutableInteger.class);
		
		MemoryIntervalTreeStorage<Transcript> transcripts = genomic.getTranscripts();
		
		IterateIntoSink<ImmutableReferenceGenomicRegion<MutableInteger>> sink = new IterateIntoSink<>(
				ei->out.fill(ei, context.getProgress()));
		
		AtomicInteger index = new AtomicInteger();
		
		EI.wrap(refs).parallelized(Math.min(8, nthreads), 1, ei->ei.map(ref->{
		
		
			TrimmedGenomicRegion reuse = new TrimmedGenomicRegion();
			IntervalTree<GenomicRegion,Void> openClusters = null;
			TreeMap<Integer,GenomicRegion> openClustersByStop = null;
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : reads.ei(ref) 
										.progress(context.getProgress(), (int)reads.size(ref), r->"Clustering reads on "+ref)
										.filter(r->r.getData().getMultiplicity(0)<=1)
										.checkOrder((a,b)->a.compareTo(b))
										.loop()) {
				

				
				ArrayList<ImmutableReferenceGenomicRegion<Transcript>> ctrans = transcripts.ei(read)
							.chain(transcripts.ei(read.toMutable().toOppositeStrand()))
							.filter(t->reuse.set(read.getRegion()).isCompatibleWith(t.getRegion())).list();
				
				if (read.getRegion().getNumParts()>1 && ctrans.isEmpty())
					continue;
				
				if (openClusters==null) {
					// first cluster
					openClusters = new IntervalTree<GenomicRegion, Void>(read.getReference());
					openClusters.put(read.getRegion(),null);
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				}
				else if (!openClusters.getReference().equals(read.getReference())) {
					// first cluster on this chromosome
					
					for (GenomicRegion clreg : openClusters.keySet()) 
						sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, new MutableInteger(index.getAndIncrement())));

					openClusters = new IntervalTree<GenomicRegion, Void>(read.getReference());
					openClusters.put(read.getRegion(),null);
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				} else {
					ArrayList<GenomicRegion> overlaps = openClusters
						.keys(read.getRegion().getStart(), read.getRegion().getStop())
						.filter(reg->intersects(reuse.set(reg),read.getRegion(),ctrans,maxdist)).list();
					
					if (overlaps.size()>0) {
						// merge clusters 
						GenomicRegion r = overlaps.get(0);
						for (int i=1; i<overlaps.size(); i++) {
							r = r.union(overlaps.get(i));
						}
						r = r.union(read.getRegion());
						
						for (GenomicRegion reg : overlaps) {
							openClusters.remove(reg);
							openClustersByStop.remove(reg.getStop());
						}
						
						openClusters.put(r, null);
						openClustersByStop.put(r.getStop(), r);

					} else {
						openClusters.put(read.getRegion(), null);
						openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
					}
					
					// prune interval tree (everything that ends left of the current start
					NavigableMap<Integer, GenomicRegion> head = openClustersByStop.headMap(read.getRegion().getStart(), false);
					for (GenomicRegion clreg : head.values()) 
						sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, new MutableInteger(index.getAndIncrement())));

					openClusters.keySet().removeAll(head.values());
					head.clear();
				}
				
			}
			
			if (openClusters!=null)
				for (GenomicRegion clreg : openClusters.keySet()) 
					sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, new MutableInteger(index.getAndIncrement())));
			return 1;
		})).drain();
		
		sink.finish();

		return null;
	}
	
	static boolean intersects(TrimmedGenomicRegion reg, GenomicRegion read,
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> ctrans, int maxdist) {
//		return reg.intersects(read.getRegion());
		int tol = maxdist-read.getTotalLength();
		
		if (read.getDistance(reg.getParent())<tol) return true;
		
		for (ImmutableReferenceGenomicRegion<Transcript> t : ctrans) {
			if (reg.isCompatibleWith(t.getRegion())){
				ArrayGenomicRegion iread = t.getRegion().induce(t.getRegion().intersect(read));
				ArrayGenomicRegion ireg = t.getRegion().induce(t.getRegion().intersect(reg));
				if (iread.getDistance(ireg)<tol) return true;
						
			}
		}
		
		return false;
	}

	
}