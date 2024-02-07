package executables;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

import gedi.util.io.text.LineOrientedFile;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Bam2Fastq {

	public static void main(String[] args) throws IOException {
		if (args.length<2) {
			usage();
			System.exit(1);
		}
		
		Progress progress = new NoProgress();
		String bam = null;
		String fastq = null;
		boolean unique = false;
		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = new ConsoleProgress(System.err);
			else if (args[i].equals("-u"))
				unique = true;
			else {
				bam = args[i++];
				if (i==args.length) {
					usage();
					System.exit(1);
				}
				fastq = args[i++];
				if (i!=args.length) {
					usage();
					System.exit(1);
				}
			}
		}
		
		
		SamReader sam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bam));
		progress.init();
		
		LineOrientedFile out = new LineOrientedFile(fastq);
		out.startWriting();
		
		
		SAMRecordIterator it = sam.iterator();
		while (it.hasNext()) {
			SAMRecord rec = it.next();
			int id = Integer.parseInt(rec.getReadName());
			
			if (!unique || rec.getIntegerAttribute("NH")==1)
				out.writef("@%d\n%s\n+\n%s\n", id, rec.getReadString(),rec.getBaseQualityString());
			
			progress.incrementProgress().setDescription(()->rec.getReadName()+" "+rec.getReferenceName()+":"+rec.getAlignmentStart());
		}

		progress.finish();
		out.finishWriting();
	}

	private static void usage() {
		System.out.println("Bam2Fastq [-u] [-p] <bam> <fastq> \n\n -p shows progress");
	}
	
}
