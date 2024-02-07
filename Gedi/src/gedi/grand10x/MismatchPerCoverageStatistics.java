package gedi.grand10x;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import gedi.core.reference.ReferenceSequence;
import gedi.grand10x.javapipeline.Grand10XDemultiplexProgram.SnpData;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineWriter;

public class MismatchPerCoverageStatistics {
	
	private static int MAX_COV = 10; 
	
	long[][][][] mm;
	
	long[][][][] pos;
	
	String cattitle;
	String[] categories;
	HashMap<String, Integer> catIndex;
	
	HashMap<String,TreeMap<Integer, SnpData>> snps = new HashMap<>();
	
	
	public MismatchPerCoverageStatistics(int readlen, String[] categories, String cattitle) {
		mm = new long[categories.length][5][5][MAX_COV+1];
		pos = new long[categories.length][5][5][readlen];
		catIndex = EI.wrap(categories).indexPosition();
		this.categories = categories;
		this.cattitle = cattitle;
	}
	
	public MismatchPerCoverageStatistics emptyCopy() {
		return new MismatchPerCoverageStatistics(pos[0][0][0].length, categories, cattitle);
	}
	
	public int getReadLen() {
		return pos[0][0][0].length;
	}

	
	public void addSnps(String r, TreeMap<Integer, SnpData> snpdata) {
		TreeMap<Integer, SnpData> s = snps.computeIfAbsent(r, x->new TreeMap<Integer, SnpData>());
		for (Integer p : snpdata.keySet())  {
			SnpData there = s.get(p);
			if (there==null)
				s.put(p, snpdata.get(p));
			else
				there.add(snpdata.get(p));
		}
	}

	public double getMmFraction(String cat, char genomic, char read, int coverage) {
		return mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][coverage]/(double)mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[genomic]][coverage];
	}
	
	public double getMmFraction(String cat, char genomic, char read) {
		long m = 0;
		long c = 0;
		for (int coverage=0; coverage<mm[0][0][0].length; coverage++) {
			m+=mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][coverage];
			c+=mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[genomic]][coverage];
		}
		return m/(double)c;
	}

	public long getCov(String cat, char genomic,int coverage) {
		return mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[genomic]][coverage];
	}
	public long getMm(String cat, char genomic, char read, int coverage) {
		return mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][coverage];
	}

	public long getMm(String cat, char genomic, char read) {
		long m = 0;
		for (int coverage=0; coverage<mm[0][0][0].length; coverage++) {
			m+=mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][coverage];
		}
		return m;
	}

	private ArrayList<String> seg = new ArrayList<String>();
	public void addSegregatingAnalysis(String locs, ReadCategory cat, double tc, double ag, long umis) {
		seg.add(String.format("%s\t%s\t%d\t%.3g\t%.3g", locs,cat.toString(),umis,tc,ag));
	}

	public void writeSeg(LineWriter sout) throws IOException {
		sout.writeLine("Location\tCategory\tUMIs\tTC\tAG");
		EI.wrap(seg).print(sout);
	}
	
	public void writeCov(LineWriter out) throws IOException {
		out.writeLine(cattitle+"\tGenomic\tRead\tReads\tMismatches\tCoverage");
		for (int cat=0; cat<pos.length; cat++) 
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					for (int r=1; r<mm[cat][i][j].length; r++)
						if (i!=j) 
							out.writef("%s\t%c\t%c\t%d\t%d\t%d\n", 
									categories[cat],
									SequenceUtils.nucleotides[i],
									SequenceUtils.nucleotides[j],
									r,
									mm[cat][i][j][r],
									mm[cat][i][i][r]);
				
	}

	public void writePos(LineWriter out) throws IOException {
		out.writeLine(cattitle+"\tGenomic\tRead\tPosition\tMismatches\tCoverage");
		for (int cat=0; cat<pos.length; cat++) 
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					for (int r=0; r<pos[cat][i][j].length; r++)
						if (i!=j)
							out.writef("%s\t%c\t%c\t%d\t%d\t%d\n", 
									categories[cat],
									SequenceUtils.nucleotides[i],
									SequenceUtils.nucleotides[j],
									r,
									pos[cat][i][j][r],
									pos[cat][i][i][r]);
				
	}
	
	
	public void writeSnps(LineWriter out, double snpconv) throws IOException {
		out.writeLine("Location\tCoverage\tMismatches\tP value");
		for (String r : snps.keySet()) {
			for (Integer p : snps.get(r).keySet()) {
				SnpData snp = snps.get(r).get(p);
				out.writef("%s:%s\t%d\t%d\t%.3g\n", r,p,snp.getCov(),snp.getMm(),snp.getPval(snpconv));
			}
		}
	}
	

	public void addReadMismatch(String cat, char genomic, char read, int pos) {
		if (pos>=this.pos[0][0][0].length)
			return;
		this.pos[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][pos]++;
	}
	
	public void addReadCoverage(String cat, char genomic, int pos) {
		if (pos>=this.pos[0][0][0].length)
			return;
		this.pos[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[genomic]][pos]++;
	}

	public void addMismatch(String cat, char genomic, char read, int covered) {
		covered = Math.min(covered, mm[0][0][0].length-1);
		mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][covered]++;
	}
	
	public void addCoverage(String cat, char genomic, int covered) {
		covered = Math.min(covered, mm[0][0][0].length-1);
		mm[catIndex.get(cat)][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[genomic]][covered]++;
	}
	
	public MismatchPerCoverageStatistics add(MismatchPerCoverageStatistics other) {
		for (int cat=0; cat<pos.length; cat++) 
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					for (int r=0; r<pos[cat][i][j].length; r++)
						pos[cat][i][j][r]+=other.pos[cat][i][j][r];
		for (int cat=0; cat<mm.length; cat++) 
			for (int i=0; i<4; i++)
				for (int j=0; j<4; j++)
					for (int r=0; r<mm[cat][i][j].length; r++)
						mm[cat][i][j][r]+=other.mm[cat][i][j][r];
		
		for (String r : other.snps.keySet())
			addSnps(r, other.snps.get(r));
		seg.addAll(other.seg);
		
		return this;
	}

	


	
}