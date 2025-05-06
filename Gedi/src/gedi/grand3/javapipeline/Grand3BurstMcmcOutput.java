package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.Grand3Utils;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.processing.SubreadProcessorMismatchBuffer;
import gedi.grand3.reads.ClippingData;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.SnpData;
import gedi.grand3.targets.TargetCollection;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.sequence.DnaSequence;

public class Grand3BurstMcmcOutput<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3BurstMcmcOutput(Grand3ParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addInput(params.targetCollection);
		addInput(params.snpFile);
		addInput(params.strandnessFile);
		addInput(params.clipFile);
		addInput(params.experimentalDesignFile);
		
		addInput(params.prefix);
		addInput(params.debug);
		
		setRunFlag(params.burstMCMC);
		
		addOutput(params.burstMCMCFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		int nthreads = getIntParameter(pind++);
		Genomic genomic = getParameter(pind++);
		GenomicRegionStorage<A> reads = getParameter(pind++);
		TargetCollection targets = getParameter(pind++);
		File snpFile = getParameter(pind++); 
		File strandnessFile = getParameter(pind++); 
		File clipFile = getParameter(pind++);
		File designFile = getParameter(pind++);
		
		String prefix = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		context.getLog().info("Processing targets");
		
		SnpData masked = new SnpData(snpFile);
		Strandness strandness = Grand3Utils.getStrandness(strandnessFile);
		ClippingData clipping = ClippingData.fromFile(clipFile);
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);

		context.getLog().info("Read "+masked.size()+" SNPs!");
		context.getLog().info("Strandness: "+strandness);
		context.getLog().info("Clipping: "+clipping);
		
		targets.checkValid();
		context.getLog().info("Using the following categories for writing burstMCMC data: "+EI.wrap(targets.getCategories(c->c.useToEstimateTargetParameters())).concat(","));
		
		ReadSource<A> source = new ReadSource<>(reads, clipping, strandness, debug);
		BarcodeHandler bcHandler = new BarcodeHandler(design.getBarcodes());
		
		LineWriter out = getOutputWriter(0);
		for (BurstMcmcResult[] res : targets.iterateRegions()
			.progress(context.getProgress(), targets.getNumRegions(), r->"Processing "+r.getData())
			.parallelized(nthreads, 5, ei->ei.map(target->{
				return processTarget(genomic, bcHandler, target, source, targets, design, masked);
			}))
			.progress(context.getProgress(), targets.getNumRegions(), r->"Finished")
			.loop()) {
			
			if (EI.wrap(res).filter(r->r.umis.size()>0).count()>0)
				for (int i=0; i<res.length; i++) {
					out.writef("%s\t%s\t%s%s\n", res[i].gene,design.getFullName(i),res[i].umis.size(),res[i].mmBuilder.toString());
				}
		}
		out.close();
		
		
		
		return null;
	}


	private static class BarcodeHandler {
		int cellLength;
		HashMap<DnaSequence,Integer> bcToIndex;
		
		public BarcodeHandler(String[] barcodes) {
			cellLength = barcodes[0].length();
			bcToIndex = EI.wrap(barcodes).map(s->new DnaSequence(s)).indexPosition();
		}
		public int getCellIndex(DnaSequence bc) {
			Integer re = bcToIndex.get(bc.subSequence(0, cellLength));
			return re==null?-1:re;
		}
		public DnaSequence getUmi(DnaSequence bc) {
			return bc.subSequence(cellLength, bc.length());
		}
	}
	private static class BurstMcmcResult  {
		private String gene;
		private StringBuilder mmBuilder = new StringBuilder();
		private HashSet<DnaSequence> umis = new HashSet<DnaSequence>();
		public BurstMcmcResult(String gene) {
			this.gene = gene;
		}
	}
	private BurstMcmcResult[] processTarget(Genomic genomic, BarcodeHandler bcHandler, ImmutableReferenceGenomicRegion<String> target, ReadSource<A> source,
			TargetCollection targets, ExperimentalDesign design, SnpData masked) {
		
		ImmutableReferenceGenomicRegion<String>  currentTargetExtended = new ImmutableReferenceGenomicRegion<>(
				target.getReference(), 
				target.getRegion().extendAll(100000, 100000).intersect(0, genomic.getLength(target.getReference().getName())),
				target.getData());
		char[] targetSeq = genomic.getSequence(currentTargetExtended).toString().toUpperCase().toCharArray();

		ReferenceSequence refInd = target.getReference().toStrandIndependent();
		
		BurstMcmcResult[] umis = new BurstMcmcResult[design.getCount()];
		for (int i=0; i<umis.length; i++)
			umis[i] = new BurstMcmcResult(target.getData());
		
		
		SubreadProcessorMismatchBuffer buff = new SubreadProcessorMismatchBuffer(1);
		for (ImmutableReferenceGenomicRegion<A> sread : source.getRawReads(target).loop()) {

			// classify read
			targets.classify(target, sread, source.getStrandness(), true, buff);
			if (buff.getMode().equals(ReadCountMode.No)) 
				continue;

			char[] readSeq;
			// obtain sequence
			if (currentTargetExtended.getRegion().contains(sread.getRegion())) {
				readSeq = SequenceUtils.extractSequence(currentTargetExtended.induce(sread.getRegion()), targetSeq);
				if (!sread.getReference().getStrand().equals(currentTargetExtended.getReference().getStrand()))
					SequenceUtils.getDnaReverseComplementInplace(readSeq);
			}
			else {
				readSeq = genomic.getSequence(sread).toString().toUpperCase().toCharArray();
			}
			
			int len = sread.getRegion().getTotalLength();
			
			AlignedReadsData data = sread.getData();
			// count mismatches per distinct
			for (int d=0; d<data.getDistinctSequences(); d++) {
				
				// count totals
				for (int i=0; i<sread.getRegion().getTotalLength(); i++) {
					int gpos = sread.map(i);
					boolean indeletion = sread.getData().mapToRead1(d, i)==-1;
					if (!indeletion && !masked.isSnp(refInd, gpos)) {
						buff.increment(0,readSeq[i]);
					}
				}
				
				for (int v=0; v<data.getVariationCount(d); v++) {
					if (data.isMismatch(d, v)) {
						int gpos = sread.map(data.getMismatchPos(d, v));
						if (!masked.isSnp(refInd, gpos) && checkSequence(data.getMismatchGenomic(d, v).charAt(0),readSeq[data.getMismatchPos(d, v)],sread,readSeq,sread.getData().getVariation(d, v))) {
							buff.increment(0, data.getMismatchGenomic(d, v).charAt(0), data.getMismatchRead(d, v).charAt(0));
						}
					}
				}
				
				if (sread.getData().getNumConditions()!=1) throw new RuntimeException("Not implemented");
				if (!(sread.getData() instanceof BarcodedAlignedReadsData)) throw new RuntimeException("Not implemented");
				
				int total = buff.getTotal(0, MetabolicLabelType._4sU);
				int conv = buff.getMismatches(0, MetabolicLabelType._4sU);
				
				BarcodedAlignedReadsData bdata = (BarcodedAlignedReadsData) sread.getData();
					
				for (DnaSequence bc : bdata.getBarcodes(d, 0)) {
					int cond = bcHandler.getCellIndex(bc);
					DnaSequence umi = bcHandler.getUmi(bc);
					umis[cond].mmBuilder.append("\t").append(conv).append("/").append(total);
					umis[cond].umis.add(umi);
				}
				
				buff.reset();
			}
			
		}
		return umis;
	}
	
	private static boolean checkSequence(char mm, char genome, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> sread, char[] readSeq, AlignedReadsVariation var) {
		if (mm=='N' && genome!='A'&&genome!='C'&&genome!='G'&&genome!='T')
			return false;
		if (mm!=genome)
			throw new RuntimeException("Sequences do not match! This is a sign that references for read mapping and Grand3 are different!\n"+sread+"\nReference sequence: "+String.valueOf(readSeq)+"\n"+var);
		return true;

	}

	
	
}
