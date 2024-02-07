package gedi.slam.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.function.Predicate;
import java.util.logging.Level;
import java.util.regex.Pattern;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.slam.GeneData;
import gedi.slam.ReadData;
import gedi.slam.SlamCollector;
import gedi.slam.SlamCollector.CounterKey;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.datastructure.dataframe.DataFrame;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializableSerializer;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.tsv.formats.Csv;
import gedi.util.math.stat.factor.Factor;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;

public class SlamCollectStatistics extends GediProgram {

	
	private HashSet<String> allowedBiotypes = new HashSet<String>(Arrays.asList("protein_coding","lincRNA","antisense","IG_LV_gene","IG_V_gene","IG_V_pseudogene","IG_D_gene","IG_J_gene","IG_J_pseudogene","IG_C_gene","IG_C_pseudogene","TR_V_gene","TR_V_pseudogene","TR_D_gene","TR_J_gene","TR_J_pseudogene","TR_C_gene","synthetic"));
	
	
	public SlamCollectStatistics(SlamParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.snpFile);
		addInput(params.locations);
		addInput(params.strandnessFile);
		addInput(params.viral);
		addInput(params.all);
		addInput(params.trim5p);
		addInput(params.trim3p);
		addInput(params.numi);
		addInput(params.mode);
		addInput(params.overlap);
		addInput(params.no4sUpattern);
		addInput(params.introns);
		addInput(params.lenientOverlap);
		addInput(params.modelall);
		addInput(params.highmem);
		
		addInput(params.prefix);
		
		addOutput(params.dataFile);
		addOutput(params.mismatchFile);
		addOutput(params.mismatchPosFile);
		addOutput(params.binomFile);
		addOutput(params.doublehitFile);
		addOutput(params.binomOverlapFile);
		addOutput(params.dataExtFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException {
		
		int nthreads = getIntParameter(0);
		Genomic genomic = getParameter(1);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(2);
		File snps = getParameter(3);
		GenomicRegionStorage<NameProvider> locations = getParameter(4);
		Strandness strandness = Strandness.valueOf(EI.lines(this.<File>getParameter(5)).first());
		String viralChr = getParameter(6);
		boolean allGenes = getBooleanParameter(7);
		int trim5p = getIntParameter(8);
		int trim3p = getIntParameter(9);
		boolean numis = getBooleanParameter(10);
		ReadCountMode mode = getParameter(11);
		ReadCountMode overlapMode = getParameter(12);
		String pat = getParameter(13);
		boolean introns = getBooleanParameter(14);
		boolean lenientOverlap = getBooleanParameter(15);
		boolean modelall = getBooleanParameter(16);
		boolean highmem = getBooleanParameter(17);
		
		if (!allGenes)
			allGenes|=genomic.getGenes().ei().map(e->e.getData()).filter(e->genomic.getGeneTable("biotype").apply(e).length()>0).count()==0;
		
		String prefix = getParameter(18);
		
		
		Pattern no4sUPattern = Pattern.compile(pat,Pattern.CASE_INSENSITIVE);
		String[] cond = reads.getMetaDataConditions();
		int[] no4sUIndices = EI.along(cond).filterInt(i->no4sUPattern.matcher(cond[i]).find()).toIntArray();
		boolean[] no4sU = new boolean[cond.length];
		for (int i : no4sUIndices)
			no4sU[i] = true;
		
		context.getLog().info("Reading SNPs...");
		DataFrame snpDf = Csv.toDataFrame(snps.getAbsolutePath());
		
		MemoryIntervalTreeStorage<Void> masked = new MemoryIntervalTreeStorage<>(Void.class);
		if (snpDf.rows()>0)
			for (Factor f : (Factor[]) snpDf.getColumn(0).getRaw()) 
				masked.add(ImmutableReferenceGenomicRegion.parse(f.name()));

		context.getLog().info("Read "+snpDf.rows()+" SNPs!");
		
		Predicate<String> geneTest = allGenes ? (g->true) : g->allowedBiotypes.contains(genomic.getGeneTable("biotype").apply(g));
		
		
		context.getLog().info("Collecting reads for each gene (Strand: "+strandness.name()+") ...");
		
		SlamCollector collector = new SlamCollector(genomic, geneTest, reads, masked, locations,strandness,trim5p, trim3p, mode, overlapMode, no4sU, introns,lenientOverlap,modelall, highmem);
		collector.ercc();
		if (viralChr!=null) {
			collector.addReferenceCounter("viral",viralChr,false);
			collector.addReferenceCounter("cellular",viralChr,true);
		}
		LineWriter numiWriter;
		if (numis) {
			context.getLog().info("Counting nUMIs!");
			numiWriter = new LineOrientedFile(prefix+".nUMIs").write();
			numiWriter.writef("Gene");
			for (int i=0; i<cond.length; i++)
				numiWriter.writef("\t%s",cond[i]);
			numiWriter.writef("\n");
			
			collector.countNumis((g,numicounts)->{
				synchronized (numiWriter) {
					try {
						if (numicounts!=null) {
							numiWriter.writef(g.getData());
							for (int i=0; i<cond.length; i++)
								numiWriter.writef("\t%s",numicounts[i]);
							numiWriter.writef("\n");
						}
					} catch (IOException e) {
						throw new RuntimeException("Cannot write nUMIs!",e);
					}
				}
			});
		} else numiWriter = null;
			
//
//		try {
//			GediCommandline cmd = new GediCommandline();
//			cmd.addParam("genomic", genomic);
//			cmd.addParam("collector", collector);
//			cmd.addParam("genetest",geneTest);
//			cmd.read();
//		} catch (ScriptException e) {}
		
		genomic.getGenes().ei()
			.progress(context.getProgress(), (int)genomic.getGenes().size(), (r)->"Processing "+r.getData())
			.filter(r->geneTest.test(r.getData()))
			.parallelized(nthreads, 5, ei->ei.unfold(collector::collect))
			.serialize(new BinarySerializableSerializer<>(GeneData.class), new PageFileWriter(getOutputFile(0).getPath()))
		.close();
		
		if (numiWriter!=null)
			numiWriter.close();
		
		EI.wrap(collector.getMismatchCategories())
			.map(cat->collector.getExtData(cat))
			.serialize(new BinarySerializableSerializer<>(GeneData.class), new PageFileWriter(getOutputFile(6).getPath()))
			.close();
	
		
		String[] conditions = reads.getMetaDataConditions();
		
		boolean[] tf = {true,false};
		boolean[] ft = {false,true};
		// write rates
		LineWriter writer = new LineOrientedFile(getOutputFile(1).getPath()).write();
		writer.write("Category\tCondition\tOrientation\tGenomic\tRead\tCoverage\tMismatches\n");
		for (String category : collector.getMismatchCategories()) {
			for (int i=0; i<conditions.length; i++) { 
				for (boolean second : ft)
				for (char g : SequenceUtils.valid_nucleotides)
					for (char r : SequenceUtils.valid_nucleotides)
						if (g!=r) {
							writer.writef("%s\t%s\t%s\t%c\t%c\t%.0f\t%.0f\n",category,conditions[i],
									second?"Second":"First",
									g,r,
									collector.getCoverage(category,i,g,second),
									collector.getMismatches(category,i,g,r,second));
						}
//				writer.writef("%s\tIntronic T\tC\t%.6g\n",conditions[i],collector.getMismatchRate("Intronic",i,'T','C'));
			}
		}

		writer.close();
		
		
		context.getLog().info("Writing statistics");
		// write rates
		int maxPos = collector.getTotalReadLength();
		
		writer = new LineOrientedFile(getOutputFile(2).getPath()).write();
		writer.write("Category\tGenomic\tRead\tPosition\tOverlap\tOpposite\tCoverage\tMismatches\n");
		for (String category : collector.getMismatchCategories()) {
//				for (CounterKey k : collector.getMismatchPositionKeys(category).loop()) {
				for (char g : SequenceUtils.valid_nucleotides)
					for (char r : SequenceUtils.valid_nucleotides) 
						if (g!=r) {
							for (boolean overlap : tf)
								for (boolean opposite : tf)
									for (int pos = 0; pos<maxPos; pos++) {
										CounterKey k = new CounterKey(SequenceUtils.inv_nucleotides[g], SequenceUtils.inv_nucleotides[r], pos, overlap, opposite);
										double cov = collector.getCoverage(category,k);
										if (cov>0) {
											double mis = collector.getMismatches(category,k);
											writer.writef("%s\t%c\t%c\t%d\t%d\t%d\t%.0f\t%.0f\n",category,k.getGenomic(),k.getRead(),
													k.getReadPos(),k.isOverlap()?1:0,k.isReadOppositeStrand()?1:0,
													cov,
													mis);
										}
									}
						}
		}
		writer.close();
		
		ArrayList<String> estimateTypes = new ArrayList<>();
		estimateTypes.add("");
		for (String type : collector.getCounterTypes())
			if (!type.endsWith("Sense") && !type.endsWith("Antisense"))
				estimateTypes.add(type);
		
		// write binom data
		LineWriter writer2 = new LineOrientedFile(getOutputFile(3).getPath()).write();
		writer2.write("n\td\tCondition\tType\tcount\n");
		for (String ty : estimateTypes)
			for (ReadData binom : collector.getBinomData(ty)) {
				AutoSparseDenseDoubleArrayCollector count = binom.getCount();
				count.process((ind,val)->{
					try {
						writer2.writef("%d\t%d\t%s\t%s\t%.0f\n", 
								binom.getTotal(),binom.getConversions(),
								conditions[ind], ty, 
								val);
					} catch (IOException e) {
						throw new RuntimeException("Cannot write statistics!",e);
					}
					return val;
				});
			}
		writer2.close();
		
		// write binom-overlap data
		LineWriter writer3 = new LineOrientedFile(getOutputFile(5).getPath()).write();
		writer3.write("n\td\tCondition\tType\tcount\n");
		for (String ty : estimateTypes)
			for (ReadData binom : collector.getBinomOverlapData(ty)) {
				AutoSparseDenseDoubleArrayCollector count = binom.getCount();
				count.process((ind,val)->{
					try {
						writer3.writef("%d\t%d\t%s\t%s\t%.0f\n", 
								binom.getTotal(),binom.getConversions(),
								conditions[ind],ty,
								val);
					} catch (IOException e) {
						throw new RuntimeException("Cannot write statistics!",e);
					}
					return val;
				});
			}
		writer3.close();
		
		
		
//		// write overlap data
//		LineWriter writer4 = new LineOrientedFile(getOutputFile(4).getPath()).write();
//		writer4.write("n\td\tCondition\tcount\n");
//		for (ReadData over : collector.getOutsideOfOverlapWithDoubleData()) {
//			AutoSparseDenseDoubleArrayCollector count = over.getCount();
//			count.process((ind,val)->{
//				try {
//					writer4.writef("%d\t%d\t%s\t%.0f\n", 
//							over.getTotal(),over.getConversions(),
//							conditions[ind], val);
//				} catch (IOException e) {
//					throw new RuntimeException("Cannot write statistics!",e);
//				}
//				return val;
//			});
//		}
//		writer4.close();
		
		// write doublehit data
		writer = new LineOrientedFile(getOutputFile(4).getPath()).write();
		writer.write("Category\tCondition\tGenomic\tRead\tCoverage\tHits\n");
		for (String category : collector.getMismatchCategories()) {
			for (int i=0; i<conditions.length; i++) { 
				for (char g : SequenceUtils.valid_nucleotides)
					for (char r : SequenceUtils.valid_nucleotides) 
						if (g!=r) {
							writer.writef("%s\t%s\t%c\t%c\t%.0f\t%.0f\n",category,conditions[i],g,r,
									collector.getDoubleHitCoverage(category,i,g),
									collector.getDoubleHits(category,i,g,r));
						}
							
//				writer.writef("%s\tIntronic T\tC\t%.6g\n",conditions[i],collector.getMismatchRate("Intronic",i,'T','C'));
			}
		}
		writer.close();
		
		
		try {
			context.getLog().info("Running R scripts for plotting");
			RRunner r = new RRunner(prefix+".plotmismatches.R");
			r.set("prefix",prefix);
			r.addSource(getClass().getResourceAsStream("/resources/R/plotmm.R"));
			r.run(true);
		} catch (Throwable e) {
			context.getLog().log(Level.SEVERE, "Could not plot!", e);
		}
		
		return null;
	}


	
	
}
