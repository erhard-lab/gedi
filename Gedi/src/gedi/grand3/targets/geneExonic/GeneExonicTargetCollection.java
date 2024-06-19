package gedi.grand3.targets.geneExonic;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.function.BiPredicate;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
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
public class GeneExonicTargetCollection implements TargetCollection {
	
	private Genomic genomic; // only for getting a target by name
	private LinkedHashMap<String,CompatibilityCategory> exonic = new LinkedHashMap<>();
	private LinkedHashMap<String,CompatibilityCategory> intronic = new LinkedHashMap<>();
	private CompatibilityCategory[] categories;

	private MemoryIntervalTreeStorage<String> genes;
	private MemoryIntervalTreeStorage<Transcript> transcripts;
	
	private Predicate<String> useForGlobal;
	private ReadCountMode mode;
	private ReadCountMode overlap;
	private boolean useExonForGlobal;
	private boolean useIntronForGlobal;
	private boolean outputExon;
	private boolean outputIntron;
	private int clipForCompatibility;
	private Function<ReferenceSequence,String> genome;
	private ToIntFunction<String> chrLength;
	private HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<Transcript>>> gene2Trans;
	private HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<Transcript>>> gene2TransSense;
	private HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<Transcript>>> gene2TransAntisense;
	private HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<Transcript>>> gene2TransBoth;
	
	public GeneExonicTargetCollection(Genomic genomic, MemoryIntervalTreeStorage<String> genes, Function<ReferenceSequence,String> genome,ToIntFunction<String> chrLength,
			Predicate<String> useForGlobal,
			boolean useExonForGlobal, boolean useIntronForGlobal, boolean outputExon, boolean outputIntron,
			MemoryIntervalTreeStorage<Transcript> transcripts, ReadCountMode mode, ReadCountMode overlap,
			int clipForCompatibility) {
		this.genomic = genomic;
		this.genes = genes;
		this.genome = genome;
		this.useForGlobal = useForGlobal;
		this.transcripts = transcripts;
		this.chrLength = chrLength;
		this.mode = mode;
		this.overlap = overlap;
		this.clipForCompatibility = clipForCompatibility;
		this.useExonForGlobal = useExonForGlobal;
		this.useIntronForGlobal = useIntronForGlobal;
		this.outputExon = outputExon;
		this.outputIntron = outputIntron;
		
		
		ArrayList<String> genomeNames = genes.ei().map(r->genome.apply(r.getReference())).unique(false).list();
		int ind=0;
		for (String g : genomeNames) {
			exonic.put(g, new CompatibilityCategory("Exonic ("+g+")", ind, true, useExonForGlobal&&useForGlobal.test(g), outputExon, false));
			intronic.put(g, new CompatibilityCategory("Intronic ("+g+")", ind+genomeNames.size(), true, useIntronForGlobal&&useForGlobal.test(g), outputIntron, false));
			ind++;
		}
		
		gene2Trans = transcripts.ei().indexMulti(t->t.getData().getGeneId());
		
		gene2TransSense = new HashMap<>();
		gene2TransAntisense = new HashMap<>();
		gene2TransBoth = new HashMap<>();
		for(ImmutableReferenceGenomicRegion<String> gene : genes.ei().loop()) {
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> s = transcripts.ei(gene).list();
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> r = transcripts.ei(gene.toOppositeStrand()).list();
			gene2TransSense.put(gene.getData(),s);
			gene2TransAntisense.put(gene.getData(),r);
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> b = new ArrayList<ImmutableReferenceGenomicRegion<Transcript>>();
			b.addAll(s); b.addAll(s);
			gene2TransBoth.put(gene.getData(), b);
		}
		
		
		categories = EI.wrap(exonic.values()).chain(EI.wrap(intronic.values())).toArray(CompatibilityCategory.class);
	}
	
	@Override
	public TargetCollection create(ReadCountMode mode, ReadCountMode overlap) {
		return new GeneExonicTargetCollection(genomic,genes, genome, chrLength, useForGlobal, useExonForGlobal, useIntronForGlobal, outputExon, outputIntron, transcripts, mode, overlap, clipForCompatibility);
	}
	
	public ImmutableReferenceGenomicRegion<String> getRegion(String name) {
		if (genomic==null)
			throw new RuntimeException("Can only get target by name if genomic object is there!");
		ImmutableReferenceGenomicRegion target = ImmutableReferenceGenomicRegion.parse(genomic, name);
		return extend1kb(target, chrLength);
	}
	
	public ExtendedIterator<ImmutableReferenceGenomicRegion<String>> iterateRegions() {
//		return genes.ei(Chromosome.obtain("C5X+")).chain(genes.ei(Chromosome.obtain("C5X-"))); 
		return genes.ei().map(g->extend1kb(g,chrLength));
	}
	
	public int getNumRegions() {
//		return iterateRegions().countInt();
		return (int) genes.size();
	}
	
	private boolean isExonic(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, ImmutableReferenceGenomicRegion<Transcript> t) {
		int inIntron = read.getRegion().intersectLengthInvert(t.getRegion());
		if (t.getData().isCoding() && t.getData().get5UtrLength(t.getReference(), t.getRegion())==0 && t.getData().get3UtrLength(t.getReference(),t.getRegion())==0) {
			// incomplete transcript
			return inIntron<=clipForCompatibility && read.getRegion().intersects(t.getRegion());
		}
		
		return inIntron<=clipForCompatibility && t.getRegion().getStart()<=read.getRegion().getStart() && t.getRegion().getEnd()>=read.getRegion().getEnd();
	}
	
	private boolean isIntronic(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, ImmutableReferenceGenomicRegion<Transcript> t) {
		int inIntron = read.getRegion().intersectLengthInvert(t.getRegion());
		if (t.getData().isCoding() && t.getData().get5UtrLength(t.getReference(), t.getRegion())==0 && t.getData().get3UtrLength(t.getReference(),t.getRegion())==0) {
			// incomplete transcript
			return inIntron>clipForCompatibility && read.getRegion().intersects(t.getRegion());
		}
		
		return inIntron>clipForCompatibility && t.getRegion().getStart()<=read.getRegion().getStart() && t.getRegion().getEnd()>=read.getRegion().getEnd();
	}
	
	private int countFilteredGenesOld(String target, Strandness strandness,
			ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, 
			boolean isStrandCorrected,
			BiPredicate<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>,ImmutableReferenceGenomicRegion<Transcript>> checker) {
		ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> cand;
		if (isStrandCorrected) 
			cand = genes.ei(read).iff(strandness.equals(Strandness.Unspecific), ei->ei.chain(genes.ei(read.toOppositeStrand()))).unfold(g->EI.wrap(gene2Trans.get(g.getData())));
		else {
			switch (strandness) {
			case Sense:cand = genes.ei(read).unfold(g->EI.wrap(gene2Trans.get(g.getData())));break;
			case Antisense: cand = genes.ei(read.toOppositeStrand()).unfold(g->EI.wrap(gene2Trans.get(g.getData()))); break;
			default: cand = genes.ei(read).chain(genes.ei(read.toOppositeStrand())).unfold(g->EI.wrap(gene2Trans.get(g.getData()))); break;
			}
		}
		cand = cand.filter(t->checker.test(read,t));
		
//		ExtendedIterator<Object> it = cand.map(t->t.getReference().getStrand().toString()+t.getData().getGeneId());
		
		// return the number of genes
		// zero if target not among them
		// negative counts if read is antisense
		int found = 0;
		HashSet<String> genes = new HashSet<>();
		while (cand.hasNext()) { 
			ImmutableReferenceGenomicRegion<Transcript> t = cand.next();
			if (t.getData().getGeneId().equals(target))
				found=t.getReference().getStrand().equals(read.getReference().getStrand())?1:-1;
			genes.add(t.getData().getGeneId());
		}
		return genes.size()*found;
	}
	
	private int countFilteredGenes(String target, Strandness strandness,
			ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, 
			boolean isStrandCorrected,
			BiPredicate<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>,ImmutableReferenceGenomicRegion<Transcript>> checker) {
		ArrayList<ImmutableReferenceGenomicRegion<Transcript>> cand;
		if (isStrandCorrected) 
			cand = strandness.equals(Strandness.Unspecific)?gene2TransBoth.get(target):gene2TransSense.get(target);
		else {
			switch (strandness) {
			case Sense:cand = gene2TransSense.get(target);break;
			case Antisense: cand = gene2TransAntisense.get(target); break;
			default: cand = gene2TransBoth.get(target); break;
			}
		}
		if (cand==null) return 0;
		
		// return the number of genes
		// zero if target not among them
		// negative counts if read is antisense
		int found = 0;
		HashSet<String> genes = new HashSet<>();
		for (ImmutableReferenceGenomicRegion<Transcript> t : cand)
			if (checker.test(read, t)) {
				if (t.getData().getGeneId().equals(target))
					found=t.getReference().getStrand().equals(read.getReference().getStrand())?1:-1;
				genes.add(t.getData().getGeneId());	
			}
		return genes.size()*found;
	}
	
	
	public void classify(ImmutableReferenceGenomicRegion<String> target, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, Strandness strandness, boolean isStrandCorrected, 
			Grand3ReadClassified classified) {
		int found = countFilteredGenes(target.getData(), strandness, read, isStrandCorrected, this::isExonic);
		CompatibilityCategory cat = null;
		if (found!=0) cat = exonic.get(genome.apply(target.getReference()));
		else {
			found = countFilteredGenes(target.getData(), strandness, read, isStrandCorrected, this::isIntronic);
			if (found!=0) cat = intronic.get(genome.apply(target.getReference()));
		}
		
		ReadCountMode mode = this.mode;
		if (found==0)
			mode = ReadCountMode.No;
		else if (overlap!=ReadCountMode.All) {
			if (overlap==ReadCountMode.CollapseAll) {
				mode = ReadCountMode.CollapseAll;
			}
			else {
				if (Math.abs(found)>1) {
					if (overlap==ReadCountMode.Unique || overlap==ReadCountMode.CollapseUnique) 
						mode = ReadCountMode.No;
					if (overlap==ReadCountMode.Divide || overlap==ReadCountMode.Weight) {
						int ufound = Math.abs(found);
						mode = this.mode.transformCounts(c->c/ufound);
					}
				} else if (overlap==ReadCountMode.CollapseUnique)
					mode = ReadCountMode.CollapseAll;
			}
		}
		
		classified.classified(target,read,cat,mode,found>=0);

	}
	
	private static <T> ImmutableReferenceGenomicRegion<T> extend1kb(ImmutableReferenceGenomicRegion<T> target, ToIntFunction<String> chrLength) {
		return  new ImmutableReferenceGenomicRegion<>(
				target.getReference(), 
				target.getRegion().extendAll(1000, 1000).intersect(0, chrLength.applyAsInt(target.getReference().getName())),
				target.getData());
	}

	@Override
	public int getNumCategories() {
		return categories.length;
	}

	@Override
	public CompatibilityCategory getCategory(int index) {
		return categories[index];
	}
	
}
