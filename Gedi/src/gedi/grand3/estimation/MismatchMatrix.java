package gedi.grand3.estimation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import gedi.grand3.experiment.ExperimentalDesign;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;

public class MismatchMatrix {

	// sampleBarcode,subread,mmtype
	protected double[][][] mmtab;
	protected ExperimentalDesign design;
	protected String[] subreads;

	public MismatchMatrix(File mmFile, ExperimentalDesign design, String[] subreads) throws IOException {
		this.subreads = subreads;
		
		ArrayList<double[][]> tab = new ArrayList<double[][]>();

		int index = 0;
		String cond = null;
		double[][] cur = null;
		for (String[] a : EI.lines(mmFile).split('\t').skip(1).loop()) {
			if (!a[0].equals(cond)) {
				cond = a[0];
				tab.add(cur=new double[subreads.length][16]);
				if (!a[0].equals(design.getSampleNameForSampleIndex(index++)))
					throw new RuntimeException("Fatal exception during reading mismatch frequency file!");
			}
			
			int s = Integer.parseInt(a[1]);
			char g = a[2].charAt(0);
			char r = a[3].charAt(0);
			double freq = Double.parseDouble(a[4]);
			cur[s][mmIndex(g,r)] = freq;
		}
		
		mmtab = tab.toArray(new double[0][][]);
		this.design = design;
	}
	
	public double getMismatchFrequencyForSample(int sample,int subread, char genomic, char read) {
		if (mmIndex(genomic, read)>=16) return 0;
		return mmtab[sample][subread][mmIndex(genomic, read)];
	}
	public double getMismatchFrequencyForCondition(int condition,int subread, char genomic, char read) {
		if (mmIndex(genomic, read)>=16) return 0;
		return mmtab[design.getIndexToSampleId()[condition]][subread][mmIndex(genomic, read)];
	}
	
	
	protected static final int mmIndex(char g, char r) {
		return SequenceUtils.inv_nucleotides[g]*4+SequenceUtils.inv_nucleotides[r];
	}
	protected static final char genomic(int mmIndex) {
		return SequenceUtils.nucleotides[mmIndex/4];
	}
	protected static final char read(int mmIndex) {
		return SequenceUtils.nucleotides[mmIndex%4];
	}

}
