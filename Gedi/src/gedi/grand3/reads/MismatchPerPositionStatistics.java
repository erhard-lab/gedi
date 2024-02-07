package gedi.grand3.reads;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.IntUnaryOperator;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.subreads.GeneSpecificUmiSenseToSubreadsConverter;
import gedi.core.data.reads.subreads.MismatchReporterWithSequence;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.Grand3Utils;
import gedi.grand3.targets.CompatibilityCategory;
import gedi.grand3.targets.TargetCollection;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;
import gedi.util.functions.IntToBooleanFunction;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.sequence.DnaSequence;

public class MismatchPerPositionStatistics extends MismatchReporterWithSequence<MismatchPerPositionStatistics> {
	
	
	private TargetCollection targets;
	private boolean debug = false;
	private int numCond;
	
	/**
	 * This is stateful and not thread safe. For parallel iteration, spawn copies per thread, and then collect them in the end!
	 * 
	 * This does not use the reindex function (to map barcodes to samples or libraries), because the reads that are reported here are the original reads!
	 * @param g
	 * @param targets 
	 */
	public MismatchPerPositionStatistics(Genomic g, TargetCollection targets, int numCond) {
		super(g);
		this.targets = targets;
		this.numCond = numCond;
		
		pCounter = new PositionCounter[catIndex(targets.getCategories()[targets.getCategories().length-1], true, true, true)+1];
	}
	
	public MismatchPerPositionStatistics setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	@Override
	public MismatchPerPositionStatistics spawn(int i) {
		return new MismatchPerPositionStatistics(genomic,targets,numCond).setDebug(debug);
	}
	@Override
	public void integrate(MismatchPerPositionStatistics other) {
		PositionCounter[] o = other.pCounter;
		for (int x=0;x<pCounter.length; x++) {
			if (pCounter[x]==null)
				pCounter[x] = o[x];
			else
				pCounter[x].add(o[x]);
		}
		this.rl1 = Math.max(rl1, other.rl1);
		this.rl2 = Math.max(rl2, other.rl2);
	}


	private ReadCountMode mode;
	private CompatibilityCategory cat;
	
	
	private class PositionCounter {
//		private CompatibilityCategory cat;
//		private boolean sense;
//		private boolean retained;
//		private boolean inOverlap;
		
		// pos, idx, cond
		// idx encodes the mismatch
		// for g=r, retained=false is empty!
		private long[][][] count;
		
		public PositionCounter(int len) {
			this.count = new long[len][16][numCond];
		}
		public void resize(int len) {
			long[][][] ncount = new long[len][][];
			for (int i=0; i<count.length; i++)
				if (i>=ncount.length)
					System.out.println();
				else
					ncount[i] = count[i];
			for (int i=count.length; i<ncount.length; i++)
				ncount[i] = new long[16][numCond];
			count = ncount;
		}
		public void add(PositionCounter o) {
			if (o==null) return;
			
			if (o.count.length>count.length)
				resize(o.count.length);
			
			for (int p=0; p<o.count.length; p++)
				for (int m=0; m<o.count[p].length; m++)
					for (int c=0; c<o.count[p][m].length; c++)
						count[p][m][c]+=o.count[p][m][c];
		}
		public void count(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, int distinct,
				ReadCountMode mode, int pos, int rl, int idx) {
			if (pos>=count.length)
				resize(pos+1); 
			
			read.getData().addCountsForDistinctInt(distinct, count[pos][idx],mode);
		}
		public long get(int p, int condition, int idx) {
			return count[p][idx][condition];
		}
		
	}
	private static final int mmIndex(char g, char r) {
		if (SequenceUtils.inv_nucleotides[g]<0 || SequenceUtils.inv_nucleotides[g]>3 || SequenceUtils.inv_nucleotides[r]<0 || SequenceUtils.inv_nucleotides[r]>3)
			return -1;
		return SequenceUtils.inv_nucleotides[g]*4+SequenceUtils.inv_nucleotides[r];
	}
	
	private static final int catIndex(CompatibilityCategory cat, boolean sense, boolean retained, boolean overlap) {
		int re = 0;
		if(overlap) re|=1<<0;
		if (retained) re|=1<<1;
		if (sense) re|=1<<2;
		re|=cat.id()<<3;
		return re;
	}
	
	private PositionCounter[] pCounter;
	private int rl1 = 0;
	private int rl2 = 0;
	
	@Override
	public void startDistinct(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read,
			int distinct, IntToBooleanFunction overlapperGetter) {
		super.startDistinct(read, distinct,overlapperGetter);
		
		boolean sense = read.getReference().getStrand().equals(currentTarget.getReference().getStrand());
		
		targets.classify(currentTarget, read, Strandness.Unspecific, false, (target,iread,cat,mode,ssense)->{
			MismatchPerPositionStatistics.this.mode = mode;
			MismatchPerPositionStatistics.this.cat = cat;
		});
		if (mode.equals(ReadCountMode.No))
			return;
		
		// update count arrays, if necessary
		int nrl1 = read.getData().hasGeometry()?read.getData().getReadLength1(distinct):read.getData().getReadLength(distinct, read.getRegion().getTotalLength());
		int nrl2 = read.getData().hasGeometry()?read.getData().getReadLength2(distinct):0;
		if (nrl1>rl1 || nrl2>rl2) {
			rl1 = Math.max(this.rl1, nrl1);
			rl2 = Math.max(this.rl2, nrl2);
			for (int i=0; i<pCounter.length; i++)
				if (pCounter[i]==null)
					pCounter[i] = new PositionCounter(rl1+rl2);
				else
					pCounter[i].resize(rl1+rl2);
		}
		
		if (debug) {
			System.out.println("Positions:");
			System.out.println("Mode: "+mode+" Sense: "+sense+" Category: "+cat);
			System.out.println(currentDistinct+" "+currentRead);
		}
		// count coverage
		for (int p=0; p<readSeq.length; p++) {
			boolean in1 = !read.getData().hasGeometry() || read.getData().isPositionInFirstRead(distinct, p);
			boolean in2 = read.getData().hasGeometry() && read.getData().isPositionInSecondRead(distinct, p);
			if (in1) {
				int rp = read.getData().mapToRead1(distinct, p);
				int idx = mmIndex(readSeq[p],readSeq[p]);
				if (idx>=0 && rp>=0) { // rp could be <0 if position is in a deletion!
					int cidx = catIndex(cat, sense, true, overlapperGetter.applyAsInt(p));
					pCounter[cidx].count(read,distinct,mode,rp,rl1+rl2,idx);
				}
			}
			// that's right, need to count twice if in overlap!
			if (in2) {
				int rp = read.getData().mapToRead2(distinct, p,rl1+rl2);
				int idx = mmIndex(SequenceUtils.getDnaComplement(readSeq[p]),SequenceUtils.getDnaComplement(readSeq[p]));
				if (idx>=0 && rp>=0) { // rp could be <0 if position is in a deletion!
					int cidx = catIndex(cat, sense, true,overlapperGetter.applyAsInt(p));
					pCounter[cidx].count(read,distinct,mode,rp,rl1+rl2,idx);
				}
			}
		}
	}
	
	@Override
	public void reportMismatch(int variation, boolean overlap, boolean retained) {
		if (!checkSequence(variation)) // unusual genomic position, do not consider at all!
			return;
		
		if (mode.equals(ReadCountMode.No))
			return;
		
		boolean sense = currentRead.getReference().getStrand().equals(currentTarget.getReference().getStrand());
		
		// count the reported mismatch
		int idx = mmIndex(currentRead.getData().getMismatchGenomic(currentDistinct, variation).charAt(0),currentRead.getData().getMismatchRead(currentDistinct, variation).charAt(0));
		if (idx>=0) {
			int p = currentRead.getData().mapToRead(
					currentDistinct, 
					currentRead.getData().getMismatchPos(currentDistinct, variation), 
					currentRead.getData().isVariationFromSecondRead(currentDistinct, variation), 
					rl1+rl2);
			if (p<0)
				throw new RuntimeException("Read mappings are corrupted: Encountered a mismatch in a deletion!\n"+currentRead.toString());
			int cidx = catIndex(cat, sense, retained, overlap);
			pCounter[cidx].count(currentRead,currentDistinct,mode,p,rl1+rl2,idx);
			
			if (debug) {
				System.out.println(currentRead.getData().getVariation(currentDistinct, variation)+" Overlap: "+overlap+" Retained: "+retained);
			}
		}
		
	}
	
	private static double nanToNull(double d) {
		if (Double.isNaN(d)) return 0;
		return d;
	}

	private static final boolean[] tf = {true,false};
	private static final boolean[] t = {true};
	private static final boolean[] f = {false};
	
	private long getCount(CompatibilityCategory cat, boolean[] sense, boolean[] retained, boolean[] overlap, int p, int condition, int idx) {
		long re = 0;
		for (boolean ssense : sense)
			for (boolean sretained : retained)
				for (boolean soverlap : overlap) {
					re+=pCounter[catIndex(cat, ssense, sretained, soverlap)].get(p, condition, idx);
				}
		return re;
	}
	
	private long getCount(CompatibilityCategory[] cat, boolean[] sense, boolean[] retained, boolean[] overlap, int p, int condition, int idx) {
		long re = 0;
		for (CompatibilityCategory c : cat)
			re+=getCount(c, sense, retained, overlap, p, condition, idx);
		return re;
	}
	
	private double getFrac(CompatibilityCategory cat, boolean[] sense, boolean[] retained, boolean[] overlap, int p, int condition, int idx, int cidx) {
		return nanToNull(getCount(cat, sense, retained, overlap, p, condition, idx)/(double)getCount(cat, sense, retained, overlap, p, condition, cidx));
	}
	
	private double getFrac(CompatibilityCategory[] cat, boolean[] sense, boolean[] retained, boolean[] overlap, int p, int condition, int idx, int cidx) {
		return nanToNull(getCount(cat, sense, retained, overlap, p, condition, idx)/(double)getCount(cat, sense, retained, overlap, p, condition, cidx));
	}
		
	
	public void write(CompatibilityCategory[] categories, File normal, Strandness strandness, String[] conditionNames) throws IOException {
				
		LineWriter sout = new LineOrientedFile(normal.getPath()).write();
		sout.writef("First read\tPosition\tGenomic\tRead\tSense\tCategory\tCorrected");
		for (String c : conditionNames) sout.writef("\t%s", c);
		sout.writeLine();
		
		
		int rl = rl1+rl2;
		
		long[] totalSense = new long[categories.length];
		long[] totalAntisense = new long[categories.length];
		for (CompatibilityCategory cat : categories) {
//			int senseCorrected = catIndex(cat, true, true, true);
//			int senseUncorrected = catIndex(cat, true, true, false);
//			int antisenseCorrected = catIndex(cat, false, true, true);
//			int antisenseUncorrected = catIndex(cat, false, true, false);
			for (char g : new char[] {'A','C','G','T'}) {
				int cidx = mmIndex(g,g);
				
				for (int c=0; c<conditionNames.length; c++) {
					totalSense[cat.id()]+=	getCount(cat, t, t, tf, 0, c, cidx);
							//pCounter[senseCorrected].get(0, c, cidx)+pCounter[senseUncorrected].get(0, c, cidx);
					totalAntisense[cat.id()]+=getCount(cat, f, t, tf, 0, c, cidx);
//							pCounter[antisenseCorrected].get(0, c, cidx)+pCounter[antisenseUncorrected].get(0, c, cidx);
				}
			}
		}
		long sumTotalSense = ArrayUtils.sum(totalSense);
		long sumTotalAntisense = ArrayUtils.sum(totalAntisense);
		double minFrac = 0;//0.001;
		
		for (CompatibilityCategory cat : categories) {
//			int senseRetainedOverlap = catIndex(cat, true, true, true);
//			int senseRemovedOverlap = catIndex(cat, true, false, true);// g=r is empty!
//			int senseRetainedOutside = catIndex(cat, true, true, false);
//			int senseRemovedOutside = catIndex(cat, true, false, false); // g=r is empty!
//			int antisenseRetainedOverlap = catIndex(cat, false, true, true);
//			int antisenseRemovedOverlap = catIndex(cat, false, false, true);
//			int antisenseRetainedOutside = catIndex(cat, false, true, false);
//			int antisenseRemovedOutside = catIndex(cat, false, false, false);
			
			for (char g : new char[] {'A','C','G','T'}) {
				int cidx = mmIndex(g,g);
				
				for (char r : new char[] {'A','C','G','T'}) {
					if (g==r) continue;
					
					int idx = mmIndex(g,r);
					for (int p=0; p<rl; p++) {
						if (totalSense[cat.id()]>minFrac*sumTotalSense && !strandness.equals(Strandness.Antisense)) {
							int up=p;
							if (EI.seq(0, conditionNames.length).mapToDouble(c->getCount(cat,t,t,t,up,c,cidx)).sum()>0) {
									sout.writef("%d\t%d\t%s\t%s\t1\t%s\t1", p<rl1?1:0,p>rl1?(rl1+rl2-p):(p+1),g,r,cat);
								for (int c=0; c<conditionNames.length; c++)
									sout.writef("\t%.4g", getFrac(cat,t,t,tf,p,c,idx,cidx));
								sout.writeLine();
							}
							sout.writef("%d\t%d\t%s\t%s\t1\t%s\t0", p<rl1?1:0,p>rl1?(rl1+rl2-p):(p+1),g,r,cat);
							for (int c=0; c<conditionNames.length; c++)
								sout.writef("\t%.4g",getFrac(cat, t, tf, tf, p, c, idx, cidx));
							sout.writeLine();
						}
						
						if (totalAntisense[cat.id()]>minFrac*sumTotalAntisense && !strandness.equals(Strandness.Sense)) {
							int up=p;
							if (EI.seq(0, conditionNames.length).mapToDouble(c->getCount(cat,f,t,t,up,c,cidx)).sum()>0) {
								sout.writef("%d\t%d\t%s\t%s\t0\t%s\t1", p<rl1?1:0, p>rl1?(rl1+rl2-p):(p+1),g,r,cat);
								for (int c=0; c<conditionNames.length; c++)
									sout.writef("\t%.4g", getFrac(cat,f,t,tf,p,c,idx,cidx));
								sout.writeLine();
							}
							sout.writef("%d\t%d\t%s\t%s\t0\t%s\t0", p<rl1?1:0,p>rl1?(rl1+rl2-p):(p+1),g,r,cat);
							for (int c=0; c<conditionNames.length; c++)
								sout.writef("\t%.4g",getFrac(cat, f, tf, tf, p, c, idx, cidx));
							sout.writeLine();
							
						}
					}
				}
			}
		}
		
		sout.close();
	}
	public double[][] getMatrix(char genomic, char read, boolean sense, CompatibilityCategory[] category, boolean firstRead,
			int[] conditions) {

		int cidx = mmIndex(genomic,genomic);
		int idx = mmIndex(genomic,read);
		
		IntUnaryOperator pp = firstRead?i->i:i->(rl1+rl2-1-i);
		double[][] re = new double[firstRead?rl1:rl2][conditions.length];

		if (conditions.length==0) return re;
		
//		int retainedOverlap = catIndex(category, sense, true, true);
//		int removedOverlap = catIndex(category, sense, false, true);// g=r is empty!
//		int retainedOutside = catIndex(category, sense, true, false);
//		int removedOutside = catIndex(category, sense, false, false); // g=r is empty!
		boolean[] ssense = {sense};
		
		for (int ci=0; ci<conditions.length; ci++) {
			int c = conditions[ci];
			for (int p=0; p<re.length; p++) {
				re[p][ci] = getFrac(category, ssense, tf, tf, pp.applyAsInt(p), c, idx, cidx);
//					(pCounter[retainedOverlap].get(pp.applyAsInt(p),c,idx)
//							+pCounter[removedOverlap].get(pp.applyAsInt(p),c,idx)
//							+pCounter[retainedOutside].get(pp.applyAsInt(p),c,idx)
//							+pCounter[removedOutside].get(pp.applyAsInt(p),c,idx)) 
//					/ (pCounter[retainedOverlap].get(pp.applyAsInt(p),c,cidx)
//							+pCounter[retainedOutside].get(pp.applyAsInt(p),c,cidx));
			}
			
		}
		
		return re;
		
	}
	public int getReadLength1() {
		return rl1;
	}

	public int getReadLength2() {
		return rl2;
	}

	
	public static void main(String[] args) throws IOException {

		Genomic g = Genomic.get("h.ens90");
		ImmutableReferenceGenomicRegion<Transcript> tr = g.getTranscripts().ei().skip(2400).filter(t->t.getRegion().getTotalLength()<1000).first();
		TargetCollection targets = Grand3Utils.getTargets(g, Strandness.Sense, ReadCountMode.Unique, ReadCountMode.Unique, 0);
		MismatchPerPositionStatistics mm = new MismatchPerPositionStatistics(g, targets, 1);
		mm.cacheRegion(new ImmutableReferenceGenomicRegion<>(tr.getReference(), tr.getRegion(),"ENSG00000228776"),0);
		
		AlignedReadsDataFactory fac1 = new AlignedReadsDataFactory(1).start();
		fac1.newDistinctSequence();
		fac1.setMultiplicity(1);
		fac1.setCount(0, 1, new DnaSequence[] {new DnaSequence("AACC")});
		fac1.addMismatch(3, 'T', 'C', false);
		fac1.addMismatch(5, 'T', 'C', false);
		ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r1 = ImmutableReferenceGenomicRegion.parse("1+:100-110",fac1.createBarcode());
		ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r2 = ImmutableReferenceGenomicRegion.parse("1+:100-110",fac1.createBarcode());
		
		r1 = tr.map(r1);
		r2 = tr.map(r2);
		
		System.out.println(tr);
		System.out.println(r1);
		System.out.println(r2);
		
		new GeneSpecificUmiSenseToSubreadsConverter(new String[] {"C"}, new String[] {"C"}, new String[] {"AA"})
			.convert("ENSG00000228776",tr.getReference(), new ArrayList<>(Arrays.asList(r1,r2)), mm)
			.drain();
		
		mm.write(targets.getCategories(), new File("normal.tsv"),Strandness.Sense, new String[] {"Test"});
		
		
	}
	

}
