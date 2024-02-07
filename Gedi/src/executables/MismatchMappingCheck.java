package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.logging.Logger;

import gedi.app.Gedi;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.mutable.MutableInteger;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;

public class MismatchMappingCheck {

	
	private static final Logger log = Logger.getLogger( MismatchMappingCheck.class.getName() );
	public static void main(String[] args) {
		try {
			start(args);
		} catch (UsageException e) {
			usage("An error occurred: "+e.getMessage(), e.additional);
			e.printStackTrace();
		} catch (Exception e) {
			System.err.println("An error occurred: "+e.getMessage());
			e.printStackTrace();
		}
	}
	

	private static void usage(String message, String additional) {
		System.err.println();
		if (message!=null){
			System.err.println(message);
			System.err.println();
		}
		System.err.println("MismatchMappingCheck -j <json> -reads <cit> [-c <condition>] -p <prefix> [-g <genomics>]");
		System.out.println();
		System.out.println("Creates fastq files with 1M reads and 0-15 T->C mismatches. The reads are sampled from genes (potentially antisense, then A->G");
		System.err.println();
		System.err.println("	Options:");
		System.err.println("	-j 				The json file to be used as a template");
		System.err.println("	-reads			Mapped reads");
		System.err.println("	-c	 			The condition to take from the reads");
		System.err.println("	-prefix 		The prefix to write files to");
		System.err.println("	-progress 		Progress status");
		System.err.println("");
		if (additional!=null) {
			System.err.println(additional);
			System.err.println("");
		}
	}
	
	
	private static class UsageException extends Exception {
		String additional = null;
		public UsageException(String msg) {
			super(msg);
		}
		public UsageException(String msg, String additional) {
			super(msg);
			this.additional = additional;
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
	private static String[] checkPair(String[] args, int index) throws UsageException {
		int p = args[index].indexOf('=');
		if (!args[index].startsWith("--") || p==-1) throw new UsageException("Not an assignment parameter (--name=val): "+args[index]);
		return new String[] {args[index].substring(2, p),args[index].substring(p+1)};
	}
	
	private static int checkIntParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new UsageException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(args[index]);
	}

	private static double checkDoubleParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isNumeric(re)) throw new UsageException("Must be a double: "+args[index-1]);
		return Double.parseDouble(args[index]);
	}
	
	
	private static void start(String[] args) throws UsageException, IOException, InterruptedException {
		Gedi.startup(true);
		
		int nthreads = Runtime.getRuntime().availableProcessors();
		String json = null;
		Genomic g = null;
		GenomicRegionStorage<AlignedReadsData> reads = null;
		Progress progress = new NoProgress();
		String condition = "0";
		String prefix = null;
		
		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-h")) {
				usage(null,null);
				System.exit(0);
			}
			else if (args[i].equals("-j")) {
				json = checkParam(args, ++i);
			}
			else if (args[i].equals("-g")) {
				ArrayList<String> gnames = new ArrayList<>();
				i = checkMultiParam(args, ++i, gnames);
				g = Genomic.get(gnames);
			}
			else if (args[i].equals("-reads")) {
				reads = Gedi.load(checkParam(args, ++i));
			}
			else if (args[i].equals("-c")) {
				condition = checkParam(args, ++i);
			}
			else if (args[i].equals("-p")) {
				prefix = checkParam(args, ++i);
			}
			else if (args[i].equals("-progress")) {
				progress = new ConsoleProgress(System.err);
			}
			else throw new UsageException("Unknown parameter: "+args[i]);
			
		}
		
		if (reads==null || json==null || prefix==null)
			throw new UsageException("Check your parameters!");
		
		
		DynamicObject param=DynamicObject.parseJson(FileUtils.readAllText(new File(json)));
		
		if (g==null)
			g = Genomic.get(param.getEntry("references").asMap().keySet());
		int cond = StringUtils.isInt(condition)?Integer.parseInt(condition):ArrayUtils.linearSearch(reads.getMetaDataConditions(),condition);
		
		Genomic ug = g;
		// count uniquely mapping reads
		int total = reads.ei()
					.progress(progress, (int)reads.size(), r->"Counting unique reads")
					.parallelized(nthreads, 1024, ei->EI.singleton(ei
							.filter(r->inGene(r,ug)!=0)
							.mapToInt(r->r.getData().getTotalCountForConditionInt(cond, ReadCountMode.Unique))
							.sum()))
					.mapToInt(ni->ni)
					.sum();
		
		System.out.println("Total number of reads:"+total);
		
		
		LinkedHashMap<String, DynamicObject> pmap = param.asMap();
		pmap.remove("adapter1");
		pmap.remove("adapter2");
		pmap.remove("adapter");
		DynamicObject[] ar = new DynamicObject[16];
		for (i=0; i<ar.length; i++)
			ar[i] = DynamicObject.fromMap("fastq",new File(prefix+".tc"+i+"_R1.fq.gz").getAbsolutePath(),"name","mm"+i);
		pmap.put("datasets", DynamicObject.from(ar));
		FileUtils.writeAllText(DynamicObject.from(pmap).toJson(), new File(prefix+".json"));
		
		double prob = 1E6/total;
		LineWriter[][] writers = new LineWriter[16][2];
		for (i=0; i<writers.length; i++) {
			writers[i][0] = new LineOrientedFile(prefix+".tc"+i+"_R1.fq.gz").write();
			writers[i][1] = new LineOrientedFile(prefix+".tc"+i+"_R2.fq.gz").write();
		}
		
		MutableInteger n = new MutableInteger();
		reads.ei()
				.progress(progress, (int)reads.size(), r->"Producing fastqs")
				.parallelized(nthreads, 1024, ()->new RandomNumbers(42),(ei,rnd)->ei
						.unfold(r->{
							int s = inGene(r,ug);
							if (s!=0) 
								return makeReads(r,ug,cond,writers.length,prob,s>0?'T':'A',s>0?'C':'G',rnd);
							return null;
						}).removeNulls())
				.forEachRemaining(a->{
					for (int m=0; m<writers.length; m++){
						if (a[m][0]!=null) {
							writers[m][0].writeLine2(mkFastq(a[m][0],n.N));
							writers[m][1].writeLine2(mkFastq(a[m][1],n.N));
						}
					}
					n.N++;
				});
		
		
		for (i=0; i<writers.length; i++) 
			for (int r=0; r<2; r++)
				writers[i][r].close();
	}


	private static String mkFastq(String seq, int n) {
		StringBuilder sb = new StringBuilder();
		sb.append("@").append(n).append("\n");
		sb.append(seq).append("\n");
		sb.append("+\n");
		for (int i=0; i<seq.length(); i++)
			sb.append("E");
		return sb.toString();
	}


	private static int inGene(ImmutableReferenceGenomicRegion<AlignedReadsData> r, Genomic g) {
		if (g.getTranscripts().ei(r).count()>0) return 1;
		if (g.getTranscripts().ei(r.toOppositeStrand()).count()>0) return -1;
		return 0;
	}

	private static ExtendedIterator<String[][]> makeReads(ImmutableReferenceGenomicRegion<AlignedReadsData> r, Genomic g, int cond, int maxtc, double prob, char genomic,char read, RandomNumbers rnd) {
		ArrayList<String[][]> re = null;
		String os = null;
		
		for (int d=0; d<r.getData().getDistinctSequences(); d++) {
			int n = r.getData().getCountInt(d, cond, ReadCountMode.Unique);
			if (n==0) continue;
			
			int k = rnd.getBinom(n, prob);
			if (k==0) continue;
			
			String[][] re1 = new String[maxtc+1][2];
			if (os==null) os = g.getSequence(r).toString();
			
			// reconstitute mismatches
			char[] s = os.toCharArray();
			
			// potential conversion positions
			IntArrayList cpos = new IntArrayList();
			for (int i=0; i<s.length; i++)
				if (s[i]==genomic)
					cpos.add(i);
			
			for (int m=0; m<re1.length; m++) {
				// introduce m conversions
				if (cpos.size()>=m) {
					cpos.shuffle(rnd);
					for (int ci=0; ci<m; ci++) 
						s[cpos.getInt(ci)] = read;
					
					re1[m][0] = String.valueOf(s, 0, r.getData().getGeometryOverlap(d)+r.getData().getGeometryBeforeOverlap(d));
					re1[m][1] = SequenceUtils.getDnaReverseComplement(String.valueOf(s, r.getData().getGeometryBeforeOverlap(d), r.getData().getGeometryOverlap(d)+r.getData().getGeometryAfterOverlap(d)));
					
					for (int ci=0; ci<m; ci++) 
						s[cpos.getInt(ci)] = genomic;
				}
			}
			if (re==null) re = new ArrayList<String[][]>();
			re.add(re1);
		}
		
		
		 
		return EI.wrap(re);
	}

}


