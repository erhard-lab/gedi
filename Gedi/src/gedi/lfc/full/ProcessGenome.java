package gedi.lfc.full;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.processing.old.CombinedGenomicRegionProcessor;
import gedi.core.processing.old.GenomicRegionProcessor;
import gedi.core.processing.old.OverlapMode;
import gedi.core.processing.old.ProcessorContext;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.lfc.LfcAlignedReadsProcessor;
import gedi.lfc.Count.StrandMode;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.datastructure.tree.redblacktree.SimpleInterval;
import gedi.util.mutable.MutableMonad;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;

import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.function.Consumer;

public class ProcessGenome {
	
	private GenomicRegionProcessor pr;
	boolean progress;
	boolean removeMultiMapping;
	OverlapMode mode;
	StrandMode strandspecific;
	HashSet<String> restrict;
	
	
	public ProcessGenome(GenomicRegionProcessor pr, boolean progress,
			boolean removeMultiMapping, OverlapMode mode,
			StrandMode strandspecific, HashSet<String> restrict) {
		this.pr = pr;
		this.progress = progress;
		this.removeMultiMapping = removeMultiMapping;
		this.mode = mode;
		this.strandspecific = strandspecific;
		this.restrict = restrict;
	}


	public void process(GenomicRegionStorage<AlignedReadsData> reads, GenomicRegionStorage<Transcript> annotation) throws IOException {
		header(pr);
		
		Progress pro = progress?new ConsoleProgress():new NoProgress();
		ProcessorContext context = new ProcessorContext();
		
		HashSet<String> refNames = new HashSet<String>(); 
		for (ReferenceSequence ref : annotation.getReferenceSequences())
			refNames.add(ref.getName());
		
		OverlapMode overlapMode = mode;
		boolean rmm = removeMultiMapping;
		GenomicRegionProcessor proc = pr;
		for (ReferenceSequence ref : reads.getReferenceSequences()) {
			
			
			pro.init();
			pro.setDescription(ref.toString());

			// mapping to genes
			HashMap<String,MutableReferenceGenomicRegion<String>> genesToRegon = new HashMap<String, MutableReferenceGenomicRegion<String>>();
			HashMap<String,IntervalTree<SimpleInterval, String>> geneToExons = new HashMap<String,IntervalTree<SimpleInterval, String>>();
			Consumer<? super ImmutableReferenceGenomicRegion<Transcript>> adder = rgr->{
				if (restrict!=null && !restrict.contains(rgr.getData().getGeneId()))
					return;

				MutableReferenceGenomicRegion<String> r = genesToRegon.get(rgr.getData().getGeneId());

				if (r==null){ 
					genesToRegon.put(rgr.getData().getGeneId(), 
						new MutableReferenceGenomicRegion<String>()
						.setReference(rgr.getReference())
						.setRegion(rgr.getRegion())
						.setData(rgr.getData().getGeneId()));
					geneToExons.put(rgr.getData().getGeneId(), new IntervalTree<>(rgr.getReference()));
				}else {
					if (!r.getReference().equals(rgr.getReference()))
						throw new RuntimeException(rgr.getData().getGeneId()+" is located on multiple chromosomes: "+r.getReference()+", "+rgr.getReference());
					r.setRegion(r.getRegion().union(rgr.getRegion()));
				}

				for (int p=0; p<rgr.getRegion().getNumParts(); p++)
					geneToExons.get(rgr.getData().getGeneId()).put(new SimpleInterval(rgr.getRegion().getStart(p), rgr.getRegion().getStop(p)),rgr.getData().getTranscriptId());

			};
			
			switch (strandspecific) {
			case No:
				annotation.iterateReferenceGenomicRegions(ref.toPlusStrand().correctName(refNames)).forEachRemaining(adder);
				annotation.iterateReferenceGenomicRegions(ref.toMinusStrand().correctName(refNames)).forEachRemaining(adder);
				break;
			case Yes:
				annotation.iterateReferenceGenomicRegions(ref.correctName(refNames)).forEachRemaining(adder);
				break;
			case Reverse:
				annotation.iterateReferenceGenomicRegions(ref.toOppositeStrand().correctName(refNames)).forEachRemaining(adder);
				break;
			}	
			
			IntervalTree<GenomicRegion,String> geneTree = new IntervalTree<GenomicRegion,String>(ref);
			for (String g : genesToRegon.keySet())
				geneTree.put(genesToRegon.get(g).getRegion(), g);

			HashMap<String,double[]> counter = new HashMap<String, double[]>();
			MutableMonad<String> uniqueGene = new MutableMonad<String>();

			reads.iterateMutableReferenceGenomicRegions(ref).forEachRemaining(read->{
				pro.incrementProgress();
				uniqueGene.Item = null;
				geneTree.iterateIntervalsIntersecting(read.getRegion(), t->true).forEachRemaining(e->{
					
					if (overlapMode.test(e.getKey(),geneToExons.get(e.getValue()),read.getRegion())) {
						double[] c = counter.computeIfAbsent(e.getValue(), g->new double[read.getData().getNumConditions()]);
						if (!rmm) {
							context.putValue(ProcessorContext.EXON_TREE, geneToExons.get(e.getValue()));
							process(proc,genesToRegon.get(e.getValue()),c,read,context);
						}
						else {
							if (uniqueGene.Item==null) uniqueGene.Item = e.getValue();
							else uniqueGene.Item = "";
						}
					}
					
					if (rmm && genesToRegon.containsKey(uniqueGene.Item)) {
						context.putValue(ProcessorContext.EXON_TREE, geneToExons.get(uniqueGene.Item));
						process(proc,genesToRegon.get(uniqueGene.Item),counter.computeIfAbsent(uniqueGene.Item, g->new double[read.getData().getNumConditions()]),read,context);
					}
				});
				
			});
			
			for (String g : counter.keySet()) 
				write(proc,counter.get(g),genesToRegon.get(g));
			
			
			pro.finish();
			
			
			
		}
		
		
	}
	
	
	private static void write(GenomicRegionProcessor pr, double[] c,
			MutableReferenceGenomicRegion<String> gene) {
		try {
			if (pr instanceof LfcAlignedReadsProcessor)
				write((LfcAlignedReadsProcessor)pr, c, gene);
			else {
				CombinedGenomicRegionProcessor cp = (CombinedGenomicRegionProcessor) pr;
				for (GenomicRegionProcessor p : cp)
					write(p,c,gene);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	private static void write(LfcAlignedReadsProcessor pr, double[] c,
			MutableReferenceGenomicRegion<String> gene) throws IOException {
		pr.setBuffer(c);
		pr.endRegion(gene, null);
	}

	private static void process(GenomicRegionProcessor pr, MutableReferenceGenomicRegion<String> gene, double[] c, MutableReferenceGenomicRegion<AlignedReadsData> read, ProcessorContext context) {
		try {
			if (pr instanceof LfcAlignedReadsProcessor)
				process((LfcAlignedReadsProcessor)pr, gene, c, read, context);
			else {
				CombinedGenomicRegionProcessor cp = (CombinedGenomicRegionProcessor) pr;
				for (GenomicRegionProcessor p : cp)
					process(p, gene, c, read, context);
			}
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	private static void process(LfcAlignedReadsProcessor pr, MutableReferenceGenomicRegion<String> gene, double[] c, MutableReferenceGenomicRegion<AlignedReadsData> read, ProcessorContext context) throws IOException {
		pr.setBuffer(c);
		pr.read(gene, read, context);
	}
	
	private static boolean header(GenomicRegionProcessor pr) throws IOException {
		if (pr instanceof LfcAlignedReadsProcessor) {
			((LfcAlignedReadsProcessor)pr).begin(null);
			return true;
		}
		else {
			CombinedGenomicRegionProcessor cp = (CombinedGenomicRegionProcessor) pr;
			for (GenomicRegionProcessor p : cp)
				if (header(p)) return true;
			return false;
		}
	}
}
