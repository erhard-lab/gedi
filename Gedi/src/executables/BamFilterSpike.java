package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import gedi.util.functions.EI;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BamFilterSpike {

	private static void usage(String msg) {
		if (msg!=null)
			System.err.println(msg);
		System.err.println("Usage: BamFilterSpike <in.bam> <spike_prefix> <out.bam>");
		System.exit(msg==null?0:1);
	}
	
	public static void main(String[] args) throws IOException {
		
		if (args[0].equals("-h")) usage(null);
		
		if (args.length!=3) usage("Wrong number of parameters!");
		
		File fin = new File(args[0]);
		String prefix = args[1];
		File fout = new File(args[2]);
		
		if (!fin.exists()) usage(args[0]+" does not exist!");
		if (fout.exists()) usage(args[2]+" already exists!");
		
		
		SamReader in = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(fin);
		
		ArrayList<SAMSequenceRecord> seqs = new ArrayList<SAMSequenceRecord>();
		for (SAMSequenceRecord s : in.getFileHeader().getSequenceDictionary().getSequences()) {
			if (!s.getSequenceName().startsWith(prefix))
				seqs.add(s);
		}

		SAMFileHeader header = in.getFileHeader().clone();
		header.getSequenceDictionary().setSequences(seqs);
		SAMFileWriter out = new SAMFileWriterFactory().setCreateIndex(false).makeBAMWriter(header, false, fout);
		
		
		int unique = 0;
		int multimapped = 0;
		int ambiguous = 0;
		int spike = 0;
		
		HashSet<Integer> ids = new HashSet<Integer>();
		for (SAMRecord[] b : EI.wrap(in.iterator()).multiplexUnsorted((a,b)->a.getReadName().equals(b.getReadName()), SAMRecord.class).loop()) {
			
			int id = Integer.parseInt(b[0].getReadName());
			if (!ids.add(id)) throw new RuntimeException("Call samtools collate first!");
			
			int pref = countPrefix(b,prefix);
			int nopref = countPrefix(b,"")-pref;
			
			if (nopref==1 && pref==0)
				unique++;
			else if (pref==0)
				multimapped++;
			else if (nopref>0 && pref>0)
				ambiguous++;
			else
				spike++;
			
			if (pref==0)
				EI.wrap(b).forEachRemaining(out::addAlignment);
		}
		
		in.close();
		out.close();
		
		System.out.println("Unique\t"+unique);
		System.out.println("Multi\t"+multimapped);
		System.out.println("Ambiguous\t"+ambiguous);
		System.out.println(prefix+"\t"+spike);
		
		
	}

	private static int countPrefix(SAMRecord[] a, String prefix) {
		int re = 0;
		for (SAMRecord r : a)
			if ((!r.getReadPairedFlag() || r.getFirstOfPairFlag() || r.getMateUnmappedFlag()) && r.getReferenceName().startsWith(prefix))
				re++;
		return re;
	}
	
}
