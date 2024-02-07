package gedi.grand3.test;

import java.io.IOException;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.subreads.MismatchReporterWithSequence;
import gedi.core.genomic.Genomic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.IntToBooleanFunction;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public class MismatchCheckerTest extends MismatchReporterWithSequence<MismatchCheckerTest> {

	

	private final static int longfrag = 0;
	private final static int shortfrag1 = 1;
	private final static int shortfrag2 = 2;
	private final static int shortfrag2c = 3;
	private final static int supershortfrag2 = 4;
	private final static int supershortfrag2c = 5;
	
	
	private int[][] mmcount = new int[302][6];
	private int[][] acount = new int[302][6];
	
	private int[][] srem1 = new int[152][];
	private int[][] sret1 = new int[152][];
	private int[][] srem2 = new int[152][];
	private int[][] sret2 = new int[152][];
	private int[][] sacount = new int[152][];
	
	public MismatchCheckerTest(Genomic g) {
		super(g);
		for (int i=0; i<srem1.length; i++) {
			srem1[i] = new int[i];
			sret1[i] = new int[i];
			srem2[i] = new int[i];
			sret2[i] = new int[i];
			sacount[i] = new int[i];
		}
	}
	
	@Override
	public MismatchCheckerTest spawn(int index) {
		return new MismatchCheckerTest(genomic);
	}

	@Override
	public void integrate(MismatchCheckerTest other) {
		for (int i=0; i<mmcount.length; i++) {
			for (int s=0; s<mmcount[i].length; s++) {
				mmcount[i][s]+=other.mmcount[i][s];
				acount[i][s]+=other.acount[i][s];
			}
		}
		for (int i=0; i<srem1.length; i++) {
			for (int s=0; s<srem1[i].length; s++) {
				srem1[i][s] += other.srem1[i][s];
				sret1[i][s] += other.sret1[i][s];
				srem2[i][s] += other.srem2[i][s];
				sret2[i][s] += other.sret2[i][s];
				sacount[i][s] += other.sacount[i][s];
			}
		}
	}
	
	private boolean hasonlyMM;
	
	@Override
	public void startDistinct(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, int distinct, IntToBooleanFunction overlapperGetter) {
		super.startDistinct(read, distinct,overlapperGetter);
		
		hasonlyMM=true;
		for (int v=0; v<read.getData().getVariationCount(distinct); v++)
			if (!read.getData().isMismatch(distinct, v))
				hasonlyMM=false;
		
		if (!hasonlyMM) return;
		
		for (int p=0; p<readSeq.length; p++) {
			boolean in1 = !read.getData().hasGeometry() || read.getData().isPositionInFirstRead(distinct, p);
			boolean in2 = read.getData().hasGeometry() && read.getData().isPositionInSecondRead(distinct, p);
			if (in1) {
				int rp = read.getData().mapToRead1(distinct, p);
				if (readSeq[p]=='A') {
					if (rp>=0) {
						int idx = longfrag;
						if (read.getData().getGeometryOverlap(distinct)>0) {
							if (read.getData().getGeometryOverlap(distinct)==read.getRegion().getTotalLength()) idx = supershortfrag2;
							else if (in1 && in2) idx = shortfrag2;
							else idx = shortfrag1;
						}
						acount[rp][idx]+=read.getData().getCount(distinct, 4, ReadCountMode.CollapseUnique);
					}
				}
			}
			if (in2) {
				int rp = read.getData().mapToRead2(distinct, p,302);
				if (readSeq[p]=='A') { // I.e. mismatch g=T and r=?!!!
					if (rp>=0) {
						if (rp<151) throw new RuntimeException(read.toString());
						int idx = longfrag;
						if (read.getData().getGeometryOverlap(distinct)>0) {
							if (read.getData().getGeometryOverlap(distinct)==read.getRegion().getTotalLength()) idx = supershortfrag2;
							else if (in1 && in2) idx = shortfrag2;
							else idx = shortfrag1;
						}
						acount[rp][idx]+=read.getData().getCount(distinct, 4, ReadCountMode.CollapseUnique);
					}
				}
			}
			
			if (read.getData().getGeometryAfterOverlap(distinct)==0 && read.getData().getGeometryBeforeOverlap(distinct)==0) {
				if (readSeq[p]=='A')
					sacount[read.getData().getGeometryOverlap(distinct)][p]+=read.getData().getCount(distinct, 4, ReadCountMode.CollapseUnique);
			}
			
		}
	}
	
	@Override
	public void reportMismatch(int variation, boolean overlap, boolean retained) {
		if (!checkSequence(variation))
			return;
		
		if (!hasonlyMM) return;
		
		int p = currentRead.getData().getMismatchPos(currentDistinct, variation);
		
		boolean in1 = !currentRead.getData().hasGeometry() || currentRead.getData().isPositionInFirstRead(currentDistinct, p);
		boolean in2 = currentRead.getData().hasGeometry() && currentRead.getData().isPositionInSecondRead(currentDistinct, p);
		
		
		if (!currentRead.getData().isVariationFromSecondRead(currentDistinct, variation)) {
			int rp = currentRead.getData().mapToRead1(currentDistinct, p);
			if (currentRead.getData().getMismatchGenomic(currentDistinct, variation).charAt(0)=='A' && currentRead.getData().getMismatchRead(currentDistinct, variation).charAt(0)=='G') {
				int idx = longfrag;
				if (currentRead.getData().getGeometryOverlap(currentDistinct)>0) {
					if (currentRead.getData().getGeometryOverlap(currentDistinct)==currentRead.getRegion().getTotalLength()) {
						idx = supershortfrag2;
						if (retained) 
							idx = supershortfrag2c;
					}
					else if (in1 && in2) {
						idx = shortfrag2;
						if (retained) 
							idx = shortfrag2c;
					}
					else idx = shortfrag1;
				}
				mmcount[rp][idx]+=currentRead.getData().getCount(currentDistinct, 4, ReadCountMode.CollapseUnique);
				
				if (currentRead.getData().getGeometryAfterOverlap(currentDistinct)==0 && currentRead.getData().getGeometryBeforeOverlap(currentDistinct)==0) {
					int[][] c = retained?sret1:srem1;
					c[currentRead.getData().getGeometryOverlap(currentDistinct)][p]+=currentRead.getData().getCount(currentDistinct, 4, ReadCountMode.CollapseUnique);
				}
			}
		}
		if (currentRead.getData().isVariationFromSecondRead(currentDistinct, variation)) {
			int rp = currentRead.getData().mapToRead2(currentDistinct, p,302);
			if (currentRead.getData().getMismatchGenomic(currentDistinct, variation).charAt(0)=='T' && currentRead.getData().getMismatchRead(currentDistinct, variation).charAt(0)=='C') {
				int idx = longfrag;
				if (currentRead.getData().getGeometryOverlap(currentDistinct)>0) {
					if (currentRead.getData().getGeometryOverlap(currentDistinct)==currentRead.getRegion().getTotalLength()) {
						idx = supershortfrag2;
						if (retained) 
							idx = supershortfrag2c;
					}
					else if (in1 && in2) {
						idx = shortfrag2;
						if (retained) 
							idx = shortfrag2c;
					}
					else idx = shortfrag1;
				}
				mmcount[rp][idx]+=currentRead.getData().getCount(currentDistinct, 4, ReadCountMode.CollapseUnique);
				
				if (currentRead.getData().getGeometryAfterOverlap(currentDistinct)==0 && currentRead.getData().getGeometryBeforeOverlap(currentDistinct)==0) {
					int[][] c = retained?sret2:srem2;
					c[currentRead.getData().getGeometryOverlap(currentDistinct)][p]+=currentRead.getData().getCount(currentDistinct, 4, ReadCountMode.CollapseUnique);
				}
			}
		}
		
		
		
	}
	
	public void write() throws IOException {
		LineWriter out = new LineOrientedFile("test.tsv").write();
		out.writef("Long.mm\tShort1.mm\tShort2.mm\tShort2.c\tsShort2.mm\tsShort2.c\tLong.t\tShort1.t\tShort2.t\tsShort2.t\n");
		for (int i=0; i<mmcount.length; i++) 
			out.writef("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					mmcount[i][0],mmcount[i][1],mmcount[i][2],mmcount[i][3],mmcount[i][4],mmcount[i][5],
					acount[i][0],acount[i][1],acount[i][2],acount[i][4]);
		out.close();
		
		
		out = new LineOrientedFile("test2.tsv").write();
		out.writef("Length\tPosition\tTotal\tRem.1\tRet.1\tRem.2\tRet.2\n");
		for (int l=0; l<sacount.length; l++) {
			for (int p=0; p<sacount[l].length; p++)
				out.writef("%d\t%d\t%d\t%d\t%d\t%d\t%d\n",l,p,
					sacount[l][p],srem1[l][p],sret1[l][p],srem2[l][p],sret2[l][p]);
		}
		out.close();
		
	}

}
