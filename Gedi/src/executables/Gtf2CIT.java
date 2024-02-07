package executables;

import java.util.ArrayList;
import java.util.HashMap;

import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.StreamLineWriter;
import gedi.util.io.text.tsv.formats.GtfFileReader;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class Gtf2CIT {
	public static void main(String[] args) {
		try {
			start(args);
		} catch (UsageException e) {
			usage("An error occurred: "+e.getMessage());
			e.printStackTrace();
		} catch (Exception e) {
			System.err.println("An error occurred: "+e.getMessage());
			e.printStackTrace();
		}
	}
	
	private static void usage(String message) {
		System.err.println();
		if (message!=null){
			System.err.println(message);
			System.err.println();
		}
		System.err.println("Gtf2CIT <Options> <gtf-file>");
		System.err.println();
		System.err.println("Options:");
		System.err.println(" -fpkm\t\t\tInclude FPKM from cufflinks");
		System.err.println(" -p\t\t\tShow progress");
		System.err.println(" -h\t\t\tShow this message");
		System.err.println();
	}
	
	
	private static class UsageException extends Exception {
		public UsageException(String msg) {
			super(msg);
		}
	}
	
	
	private static int checkMultiParam(String[] args, int index, ArrayList<String> re) throws UsageException {
		while (index<args.length && !args[index].startsWith("-")) 
			re.add(args[index++]);
		return index-1;
	}
	private static String checkParam(String[] args, int index) throws UsageException {
		if (index>=args.length || args[index].startsWith("-")) throw new UsageException("Missing argument for "+args[index-1]);
		return args[index];
	}
	private static int checkIntParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new UsageException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(args[index]);
	}
	private static double checkDoubleParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isNumeric(re)) throw new UsageException("Must be an double: "+args[index-1]);
		return Double.parseDouble(args[index]);
	}
	
	public static void start(String[] args) throws Exception {
		Gedi.startup(true);
		
		boolean progress = true;
		boolean fpkm = false;
		double minfrac = 0.05;
		
		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-p")) 
				progress = true;
			else if (args[i].equals("-fpkm")) 
				fpkm = true;
			else if (args[i].equals("-minfrac")) 
				minfrac = checkDoubleParam(args, ++i);
			else if (args[i].equals("-h")) { 
				usage(null);
				System.exit(0);
			}
			else 
				break;
		}
		
		if (i+1!=args.length)
			usage("Specify gtf file!");
		
		
		ConsoleProgress prg = new ConsoleProgress(System.err);

		GtfFileReader gtf = new GtfFileReader(args[i], "exon", "CDS");
		if (progress)
			gtf.setProgress(prg);
		gtf.setTableOutput(FileUtils.getExtensionSibling(args[i],"genes.tsv"),FileUtils.getExtensionSibling(args[i],"transcripts.tsv"),"FPKM","frac");
		
		MemoryIntervalTreeStorage<Transcript> mem = gtf.readIntoMemoryTakeFirst(new StreamLineWriter(System.err));
		prg.finish();
		
		if (fpkm) {
			HeaderLine h = new HeaderLine();
			HashMap<String, Double> ftab = EI.lines(FileUtils.getExtensionSibling(args[i],"transcripts.tsv"))
					.header(h)
					.split("\t")
					.indexOverwrite(a->a[h.get("Transcript ID")],a->Double.parseDouble(a[h.get("FPKM")]));

			HashMap<String, Double> fractab = EI.lines(FileUtils.getExtensionSibling(args[i],"transcripts.tsv"))
					.header(h)
					.split("\t")
					.indexOverwrite(a->a[h.get("Transcript ID")],a->Double.parseDouble(a[h.get("frac")]));

			
			double uminfrac = minfrac;
			CenteredDiskIntervalTreeStorage<ScoreNameAnnotation> cit = new CenteredDiskIntervalTreeStorage<>(FileUtils.getExtensionSibling(args[i],"cit"),ScoreNameAnnotation.class);
			cit.fill(mem.ei()
					.filter(r->fractab.get(r.getData().getTranscriptId())>uminfrac)
					.map(r->new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), new ScoreNameAnnotation(r.getData().getTranscriptId(), ftab.get(r.getData().getTranscriptId()))))
					,progress?new ConsoleProgress(System.err):null);
		}

		
	}


	
}
