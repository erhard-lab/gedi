package gedi.grand3.targets.geneExonic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import clojure.main;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.Grand3ReadClassified;
import gedi.grand3.targets.TargetCollection;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

/**
 * Only genes from the first genomic without mito are used for global parameter estimation!
 * @author erhard
 *
 */
public class SpliceAcceptorTargetCollection extends GeneExonicTargetCollection {
	
	private int minOverlap;
	
	public SpliceAcceptorTargetCollection(Genomic genomic, MemoryIntervalTreeStorage<String> genes, Function<ReferenceSequence,String> genome,ToIntFunction<String> chrLength,
			Predicate<String> useForGlobal,
			MemoryIntervalTreeStorage<Transcript> transcripts, ReadCountMode mode, ReadCountMode overlap,
			int clipForCompatibility, int minOverlap) {
		super(genomic, genes, genome, chrLength, useForGlobal, false, true, true, true, transcripts, mode, overlap, clipForCompatibility);
		this.minOverlap = minOverlap;
	}
	
	@Override
	public TargetCollection create(ReadCountMode mode, ReadCountMode overlap) {
		return new SpliceAcceptorTargetCollection(genomic,genes, genome, chrLength, useForGlobal, transcripts, mode, overlap, clipForCompatibility,minOverlap);
	}
	
	private static class MyClassify implements Grand3ReadClassified {

		private Collection<String> targets;
		private ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read;
		private CompatibilityCategory cat;
		private ReadCountMode mode;
		private boolean sense;
		
		@Override
		public void classified(Collection<String> targets, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read,
				CompatibilityCategory cat, ReadCountMode mode, boolean sense) {
			this.targets = targets;
			this.read = read;
			this.cat = cat;
			this.mode = mode;
			this.sense = sense;
		}

		@Override
		public String toString() {
			return "MyClassify [targets=" + targets + ", read=" + read + ", cat=" + cat + ", mode=" + mode + ", sense="
					+ sense + "]";
		}
		
		
	}
	
	
	public void classify(ImmutableReferenceGenomicRegion<String> target, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, Strandness strandness, boolean isStrandCorrected, 
			Grand3ReadClassified classified) {
		MyClassify cl = new MyClassify();
		super.classify(target, read, strandness, isStrandCorrected, cl);
		
		ArrayList<ImmutableReferenceGenomicRegion<Transcript>> cand;
		if (isStrandCorrected) 
			cand = strandness.equals(Strandness.Unspecific)?gene2TransBoth.get(target.getData()):gene2TransSense.get(target.getData());
		else {
			switch (strandness) {
			case Sense:cand = gene2TransSense.get(target.getData());break;
			case Antisense: cand = gene2TransAntisense.get(target.getData()); break;
			default: cand = gene2TransBoth.get(target.getData()); break;
			}
		}
		if (cand==null) return;
		

		HashSet<String> sites = new HashSet<>();
		for (ImmutableReferenceGenomicRegion<Transcript> t : cand) {
			if (t.getReference().isPlus()) {
				int donor = t.getRegion().getStop(0);
				for (int ex=1; ex<t.getRegion().getNumParts(); ex++) {
					int acceptor = t.getRegion().getStart(ex);
					int acceptorReadPart = read.getRegion().getEnclosingPartIndex(acceptor);
					
					if (acceptorReadPart>=0) {
						int readPos = read.getRegion().induce(acceptor);
						if (readPos>=minOverlap && read.getRegion().getTotalLength()-readPos>=minOverlap) {
							boolean intronexon = read.getRegion().contains(acceptor-minOverlap);
							boolean exonexon = !intronexon && read.getRegion().getEnclosingPartIndex(donor)==acceptorReadPart-1;
							
							if (intronexon || exonexon) 
								sites.add(String.format("%s_%s%s:%d_%s", 
											target.getData(),
											read.getReference().getName(),
											read.getReference().getStrand().toString(),
											acceptor,
											intronexon?"IE":"EE"
												));
							
						}
					}
					
					donor = t.getRegion().getStop(ex);
				}
			} else {
				// negative strand!
				int donor = t.getRegion().getStart(t.getRegion().getNumParts()-1);
				for (int ex=t.getRegion().getNumParts()-2; ex>=0; ex--) {
					int acceptor = t.getRegion().getStop(ex);
					int acceptorReadPart = read.getRegion().getEnclosingPartIndex(acceptor);
					
					if (acceptorReadPart>=0) {
						int readPos = read.induce(acceptor);
						if (readPos>=minOverlap && read.getRegion().getTotalLength()-readPos>=minOverlap) {
							boolean intronexon = read.getRegion().contains(acceptor+minOverlap);
							boolean exonexon = !intronexon && read.getRegion().getEnclosingPartIndex(donor)==acceptorReadPart+1;
							
							if (intronexon || exonexon) 
								sites.add(String.format("%s_%s%s:%d_%s", 
											target.getData(),
											read.getReference().getName(),
											read.getReference().getStrand().toString(),
											acceptor,
											intronexon?"IE":"EE"
												));
							
						}
					}
					
					donor = t.getRegion().getStart(ex);
				}
				
			}
		}
		
		if (sites.size()>0) 
			classified.classified(sites, cl.read, exonic.get(genome.apply(target.getReference())), cl.mode, cl.sense);
		else
			classified.classified(cl.targets, cl.read, cl.cat, cl.mode, cl.sense);
		
	}
	
	public static void main(String[] args) {
		Genomic g = Genomic.get("h.ens90");
		String first = g.getOriginList().get(0);
		TargetCollection targets =  new SpliceAcceptorTargetCollection(
					g,
					g.getGenes(),
					r->r.isMitochondrial()?"Mito":g.getOrigin(r).getId(),
					r->g.getLength(r),
					s->s.equals(first),
					g.getTranscripts(),
					ReadCountMode.Unique,
					ReadCountMode.Unique,
					5,10);
		
		MyClassify my = new MyClassify();
		
		ImmutableReferenceGenomicRegion<String> gene = g.getGeneMapping().apply(g.getGeneTable("symbol","geneId").apply("IFITM3")).toImmutable();
		targets.classify(gene, ImmutableReferenceGenomicRegion.parse("11-:319745-320130"), Strandness.Sense, true, my);
		System.out.println(my);
	}
	
	
}
