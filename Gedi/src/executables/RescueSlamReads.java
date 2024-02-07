package executables;

import java.io.File;
import java.io.IOException;
import java.util.Comparator;

import gedi.bam.tools.BamUtils;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.region.bam.FactoryGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.EI;
import gedi.util.mutable.MutableInteger;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class RescueSlamReads {

	public static void main(String[] args) throws IOException {
		
		SamReader sam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(args[0]));
		SAMFileWriter out = new SAMFileWriterFactory().makeBAMWriter(sam.getFileHeader(), false, FileUtils.getExtensionSibling(new File(args[0]), ".rescued.bam"));
		
		MutableInteger uniques = new MutableInteger();
		MutableInteger multi = new MutableInteger();
		EI.wrap(sam.iterator())
			.multiplexUnsorted((a, b)->a.getReadName().equals(b.getReadName()),SAMRecord.class)
//			.parallelized(12, 1024, ei->ei.map(a->rescue(a)))
			.map(a->rescue(a))
			.sideEffect(a->{ if (a.length==1) uniques.N++; else if (a.length>1) multi.N++;})
			.unfold(a->EI.wrap(a))
			.forEachRemaining(s->out.addAlignment(s));
		
		System.out.println("Unique reads: "+uniques.N);
		System.out.println("Multireads: "+multi.N);
		out.close();
		sam.close();
	}
	
	private static int compBoth(int[] a1, int[] a2) {
		if (a1[1]<a2[1]) return -1;
		if (a1[1]>a2[1]) return 1;
		if (a1[0]<a2[0]) return 1;
		if (a1[0]>a2[0]) return -1;
		return 0;
	}
	
//	private static int compSum(int[] a1, int[] a2) {
//		return Integer.compare(a1[0]+a1[1], a2[0]+a2[1]);
//	}
	
	private static final int[] COND = {1};
	public static SAMRecord[] rescue(SAMRecord[] a) {
		if (a.length==1) return a;
		
		IntArrayList bests = prioritize(a,RescueSlamReads::compBoth);
//		IntArrayList bestsnotc = prioritize(a,RescueSlamReads::compSum);
//		if (bests.size()!=bestsnotc.size()) {
//			System.out.println("No T-C aware: "+bestsnotc.size()+" Aware: "+bests.size()+" "+a[0].getReadName());
//		}
		
		a = ArrayUtils.restrict(a, bests.toIntArray());
		for (int i=0; i<a.length; i++) {
			a[i].setNotPrimaryAlignmentFlag(i!=0);
			if (a[i].getIntegerAttribute("NH")!=null)
				a[i].setAttribute("NH", a.length);
		}
		
		if (bests.size()==1) 
			a[0].setMappingQuality(255);
		
		
		
		return a;
	}
	
	private static IntArrayList prioritize(SAMRecord[] a, Comparator<int[]> comp) {
		IntArrayList bests = new IntArrayList(a.length);
		int[] best = {0,9999};
		int[] current = {0,0};
		for (int i=0; i<a.length; i++) {
			FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(a[i],COND,false,false,null);
			fac.add(a[i], 0);
			DefaultAlignedReadsData ard = fac.create();
			current[0] = current[1] = 0;
			for (int v=0; v<ard.getVariationCount(0); v++) {
				if (ard.isMismatch(0, v) && ard.getMismatchGenomic(0, v).charAt(0)=='T' && ard.getMismatchRead(0, v).charAt(0)=='C')
					current[0]++;
				else
					current[1]++;
			}
			if (comp.compare(current, best)<=0) {
				if (comp.compare(current, best)<0) { 
					bests.clear();
					best[0] = current[0]; best[1] = current[1];
				}
				bests.add(i);
			}
		}
		return bests;
	}
	
}
