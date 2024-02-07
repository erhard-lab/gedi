package executables;

import java.util.ArrayList;

import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class SplitExonicCIT {
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
		System.err.println("SplitExonicCIT <Options> <cit-file>");
		System.err.println();
		System.err.println("Options:");
		System.err.println(" -g <genome1 genome2 ...>\t\t\tGenome names");
		System.err.println(" -cit <cit file>\t\t\tTranscript CITs");
		System.err.println(" -p\t\t\tShow progress");
		System.err.println(" -nthreads <n>\t\t\tNumber of threads");
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
	
	public static void start(String[] args) throws Exception {
		Gedi.startup(true);
		
		MemoryIntervalTreeStorage<Transcript> trans = new MemoryIntervalTreeStorage<Transcript>(Transcript.class);
		boolean progress = true;
		boolean nosrna = false;
		int nthreads = Runtime.getRuntime().availableProcessors();
		boolean simplify = false;
		Strandness strandness = Strandness.Unspecific;
		
		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-g")) {
				ArrayList<String> gnames = new ArrayList<>();
				i = checkMultiParam(args, ++i, gnames);
				trans.fill(Genomic.get(gnames).getTranscripts());
			}
			if (args[i].equals("-cit")) {
				trans.fill(new CenteredDiskIntervalTreeStorage<Transcript>(checkParam(args, ++i)));
			}
			else if (args[i].equals("-p")) 
				progress = true;
			else if (args[i].equals("-nosRNA")) 
				nosrna = true;
			else if (args[i].equals("-strandness")) 
				strandness = Strandness.valueOf(checkParam(args, ++i));
			else if (args[i].equals("-s")) 
				simplify = true;
			else if (args[i].equals("-nthreads")) 
				nthreads = checkIntParam(args, ++i);
			else if (args[i].equals("-h")) { 
				usage(null);
				System.exit(0);
			}
			else 
				break;
		}
		
		if (i+1!=args.length)
			usage("Specify CIT file!");
		
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> in = new CenteredDiskIntervalTreeStorage<>(args[i]); 
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> exonic = new CenteredDiskIntervalTreeStorage<>(FileUtils.getExtensionSibling(args[i],"exonic.cit"),in.getType()); 
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> other = new CenteredDiskIntervalTreeStorage<>(FileUtils.getExtensionSibling(args[i],"nonexonic.cit"),in.getType()); 
		
		IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> exonicSink = new IterateIntoSink<>(exonic::fill);  
		IterateIntoSink<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> otherSink = new IterateIntoSink<>(other::fill);  
		
		
		boolean usimplify = simplify;
		boolean unosrna = nosrna;
		Strandness ustrandness = strandness;
		in.ei()
			.iff(progress, ei->ei.progress(new ConsoleProgress(System.err), (int)in.size(), l->l.toLocationString()))
			.parallelized(nthreads, 4096, ei->ei.map(r->{
				if (usimplify)
					r = simplify(r);
				if (r==null) 
					return 1;
				if (testExonic(r,trans,ustrandness,unosrna))
					exonicSink.accept(r);
				else
					otherSink.accept(r);
				return 1;
			}))
			.drain();
		
		exonicSink.finish();
		otherSink.finish();
		
		exonic.setMetaData(in.getMetaData());
		other.setMetaData(in.getMetaData());
	}

	private static ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> simplify(
			ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> r) {

		if (r.getData().hasGeometry())
			throw new RuntimeException("No simplify for paired end!");
		
		int[] c = r.getData().getTotalCountsForConditionsInt(ReadCountMode.Unique);
		if (ArrayUtils.sum(c)==0) return null;
		return new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), AlignedReadsDataFactory.createSimple(c, r.getData().hasNonzeroInformation()));
	}

	private static boolean testExonic(ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> r,
			MemoryIntervalTreeStorage<Transcript> trans, Strandness strandness, boolean nosrna) {
		
		if (r.getData().hasGeometry())
			throw new RuntimeException("Not for for paired end!");
		if (r.getRegion().getNumParts()>1) return true;
		
		ExtendedIterator<ImmutableReferenceGenomicRegion<Transcript>> it = EI.empty();
		if (!strandness.equals(Strandness.Antisense))
			it = it.chain(trans.ei(r));
		if (!strandness.equals(Strandness.Sense))
			it = it.chain(trans.ei(r.toOppositeStrand()));
			
		if (nosrna)
			it = it.filter(t->!issRNA(t));
		
		for (ImmutableReferenceGenomicRegion<Transcript> t : it.loop()) {
			if (r.getRegion().intersect(t.getRegion()).getTotalLength()>=r.getRegion().getTotalLength()-10)
				return true;
		}
		return false;
	}
	
	public static boolean issRNA(ImmutableReferenceGenomicRegion<Transcript> t) {
		return t.getData().issRNA(t);
	}

	
}
