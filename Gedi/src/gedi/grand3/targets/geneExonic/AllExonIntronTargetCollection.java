package gedi.grand3.targets.geneExonic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
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
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPart;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.Grand3ReadClassified;
import gedi.grand3.targets.TargetCollection;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

/**
 * Only genes from the first genomic without mito are used for global parameter estimation!
 * @author erhard
 *
 */
public class AllExonIntronTargetCollection extends GeneExonicTargetCollection {
	
	
	public AllExonIntronTargetCollection(Genomic genomic, MemoryIntervalTreeStorage<String> genes, Function<ReferenceSequence,String> genome,ToIntFunction<String> chrLength,
			Predicate<String> useForGlobal,
			boolean useExonForGlobal, boolean useIntronForGlobal, boolean outputExon, boolean outputIntron,
			MemoryIntervalTreeStorage<Transcript> transcripts, ReadCountMode mode, ReadCountMode overlap,
			int clipForCompatibility) {
		super(genomic, genes, genome, chrLength, useForGlobal, useExonForGlobal, useIntronForGlobal, outputExon, outputIntron, transcripts, mode, overlap, clipForCompatibility);
		genomic.getUnionTranscriptMapping();
	}
	
	@Override
	public TargetCollection create(ReadCountMode mode, ReadCountMode overlap) {
		return new AllExonIntronTargetCollection(genomic,genes, genome, chrLength, useForGlobal,useExonForGlobal, useIntronForGlobal, outputExon, outputIntron, transcripts, mode, overlap, clipForCompatibility);
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

		HashSet<String> sites = new HashSet<>();
		ReferenceGenomicRegion<String> t = genomic.getUnionTranscriptMapping().apply(target.getData());
		ReferenceGenomicRegion<String> body = new ImmutableReferenceGenomicRegion<String>(t.getReference(), t.getRegion().removeIntrons(),t.getData());
		
		if (body.contains(read) && cl.cat!=null) {
		
			ArrayGenomicRegion rind = body.induce(read.getRegion());
			ArrayGenomicRegion ex = body.induce(t.getRegion());
			ArrayGenomicRegion intr = ex.invert();
			
			if (intr.intersects(rind)&& !ex.intersects(rind)) {
				// in intron
				int index = intr.getEnclosingPartIndex(rind.getStart());
				GenomicRegion intron = body.map(intr.getPart(index).asRegion());
				
				if (cl.cat.getName().startsWith(INTRONIC))  
					sites.add(String.format("%s_%s_%d_%s%s:%d-%d", 
						target.getData(),
						"intron",
						index,
						read.getReference().getName(),
						read.getReference().getStrand().toString(),
						intron.getStart(),intron.getEnd()
							));
			} else if (ex.intersects(rind) && !intr.intersects(rind)) {
				// in exon
				int index = ex.getEnclosingPartIndex(rind.getStart());
				GenomicRegion exon = body.map(ex.getPart(index).asRegion());
				
				
				// there are strange edge cases (e.g. two transcripts with exons like this:
				//         ###############
				//                    #################
				//              rrrrrrrrrrrrrrrrr
				// the read (rrr...) is by the cat definition intronic (as there is no transcript in which it is exonic)
				// w.r.t. the union transcript, it is exonic
				
				if (cl.cat.getName().startsWith(EXONIC))  
					sites.add(String.format("%s_%s_%d_%s%s:%d-%d", 
							target.getData(),
							"exon",
							index,
							read.getReference().getName(),
							read.getReference().getStrand().toString(),
							exon.getStart(),exon.getEnd()
								));
			}
		}
		
		classified.classified(sites, cl.read, cl.cat, cl.mode, cl.sense);
		
	}
	
	public static void main(String[] args) {
		Genomic g = Genomic.get("h.ens90");
		String first = g.getOriginList().get(0);
		TargetCollection targets =  new AllExonIntronTargetCollection(
					g,
					g.getGenes(),
					r->r.isMitochondrial()?"Mito":g.getOrigin(r).getId(),
					r->g.getLength(r),
					s->s.equals(first),
					false,true,true,true,
					g.getTranscripts(),
					ReadCountMode.Unique,
					ReadCountMode.Unique,
					5);
		
		MyClassify my = new MyClassify();
		
		ImmutableReferenceGenomicRegion<String> gene = g.getGeneMapping().apply(g.getGeneTable("symbol","geneId").apply("IFITM3")).toImmutable();
		targets.classify(gene, ImmutableReferenceGenomicRegion.parse("11-:319745-320130"), Strandness.Sense, true, my);
		System.out.println(my);
		targets.classify(gene, ImmutableReferenceGenomicRegion.parse("11-:319745-319980"), Strandness.Sense, true, my);
		System.out.println(my);
		targets.classify(gene, ImmutableReferenceGenomicRegion.parse("11-:320060-320070"), Strandness.Sense, true, my);
		System.out.println(my);
	}
	
	
}
