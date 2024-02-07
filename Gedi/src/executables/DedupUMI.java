package executables;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.function.BiFunction;
import java.util.logging.Level;
import java.util.logging.Logger;

import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsDeletion;
import gedi.core.data.reads.AlignedReadsMismatch;
import gedi.core.data.reads.AlignedReadsSoftclip;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.feature.output.PlotReport;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.graph.SimpleDirectedGraph;
import gedi.util.datastructure.graph.SimpleDirectedGraph.AdjacencyNode;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.counting.Counter;
import gedi.util.mutable.MutableInteger;
import gedi.util.r.RRunner;
import gedi.util.sequence.DnaSequence;
import gedi.util.sequence.MismatchGraphBuilder;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class DedupUMI {

	private static final Logger log = Logger.getLogger( DisplayRMQ.class.getName() );
	public static void main(String[] args) {
		try {
			start(args);
		} catch (UsageException e) {
			usage("An error occurred: "+e.getMessage());
			if (ArrayUtils.find(args, "-D")>=0)
				e.printStackTrace();
		} catch (Exception e) {
			System.err.println("An error occurred: "+e.getMessage());
			if (ArrayUtils.find(args, "-D")>=0)
				e.printStackTrace();
		}
	}
	
	private static void usage(String message) {
		System.err.println();
		if (message!=null){
			System.err.println(message);
			System.err.println();
		}
		System.err.println("DedupUMI <Options> file.cit");
		System.err.println();
		System.err.println("Options:");
		System.err.println(" -prefix <prefix>\t\t\tprefix for .cit, .dedup.tsv, dedupmm.tsv");
		System.err.println(" -min <n>\t\t\tMinimal number of reads per UMI (use 2 if you want to discard UMIs with only a single read) default=1");
		System.err.println(" -stats \t\t\tProduce stats");
		System.err.println(" -donoumi \t\t\tForce demultiplexing (of all reads mapping to a location) without UMIs");
		System.err.println(" -plot \t\t\tProduce stats and plots");
		System.err.println(" -nocollapse \t\t\tDo not collapse reads (i.e. no deduplication, only remove the barcodes)");
		System.err.println(" -mm \t\t\tDo the same mismatch correction as UMItools");
		System.err.println();
		System.err.println(" -h\t\t\tShow this message");
		System.err.println(" -progress\t\t\tShow progress");
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
		if (!StringUtils.isNumeric(re)) throw new UsageException("Must be a double: "+args[index-1]);
		return Double.parseDouble(args[index]);
	}
	
	public static void start(String[] args) throws Exception {
		Gedi.startup(false);
		
		boolean donoumi = false;
		boolean nocollapse = false;
		boolean stats = false;
		boolean plot = false;
		boolean progress = false;
		boolean mismatchCorrection = false;
		String prefix = null;
		int minReadCount = 1;
		
		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-h")) {
				usage(null);
				return;
			}
			else if (args[i].equals("-prefix")) 
				prefix = checkParam(args,++i);
			else if (args[i].equals("-stats")) 
				stats = true;
			else if (args[i].equals("-nocollapse"))
				nocollapse=true;
			else if (args[i].equals("-mm"))
				mismatchCorrection=true;
			else if (args[i].equals("-min"))
				minReadCount=checkIntParam(args, ++i);
			else if (args[i].equals("-noumi"))
				donoumi=true;
			else if (args[i].equals("-plot")) 
				stats = plot = true;
			else if (args[i].equals("-progress")) 
				progress = true;
			else if (args[i].equals("-D")) {} 
			else if (!args[i].startsWith("-")) 
					break;
			else throw new UsageException("Unknown parameter: "+args[i]);
		}
		
		String inp = checkParam(args, i++);

		if (!new File(inp).exists()) {
			usage("File "+inp+" does not exist!");
			System.exit(1);
		}

		CenteredDiskIntervalTreeStorage<? extends AlignedReadsData> storage = new CenteredDiskIntervalTreeStorage<>(inp);
		
		boolean replace = inp.equals(prefix+".cit");
		String outfile = replace?prefix+".tmp.cit":prefix+".cit";
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData>(outfile,DefaultAlignedReadsData.class,false);

		int[][][] mm = new int[5][5][3];
		IntArrayList dc = new IntArrayList(100);
		
		ConsoleProgress prog = progress?new ConsoleProgress(System.err):null;
		
		boolean unocollapse = nocollapse;
		boolean nmismatchCorrection = mismatchCorrection;
		int uminReadCount = minReadCount;
		
		if (storage.getRandomRecord().getClass()==BarcodedAlignedReadsData.class) {
			out.fill(((CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData>)storage).ei()
					.iff(progress, ei->ei.progress(prog, (int)storage.size(), l->l.toLocationString()))
					.map(r->dedupUmi(r,dc,mm,unocollapse, nmismatchCorrection,uminReadCount))
					.removeNulls()
					,prog
					);
		} else if (donoumi) {
			out.fill(storage.ei()
					.iff(progress, ei->ei.progress(prog, (int)storage.size(), l->l.toLocationString()))
					.map(r->dedupNoUmi(r,dc,mm,unocollapse, nmismatchCorrection,uminReadCount))
					.removeNulls()
					,prog
					);
		} else {
			throw new UsageException("No UMIs found in file; if you really want to do this, start with -noumi!");
		}
		
		
		if (replace) {
			out=null;
			new File(outfile).renameTo(new File(prefix+".cit"));
		} else {
			DynamicObject meta = storage.getMetaData();
			out.setMetaData(meta);
			out=null;
		}
		
		if (stats) {
			LineWriter dedup = new LineOrientedFile(prefix+".dedup.count.tsv").write();
			dedup.writeLine("Duplication\tCount");
			for (int d=1; d<dc.size(); d++)
				if (dc.getInt(d)>0)
					dedup.writef("%d\t%d\n", d,dc.getInt(d));
			dedup.close();
			
			LineWriter dedupmm = new LineOrientedFile(prefix+".dedup.mm.tsv").write();
			dedupmm.writeLine("Genomic\tRead\tRetainedDedup\tRetainedDup\tTotal");
			for (int g=0; g<mm.length; g++) for (int r=0; r<mm[g].length; r++)
				if (mm[g][r][2]>0)
					dedupmm.writef("%s\t%s\t%d\t%d\t%d\n",SequenceUtils.nucleotides[g],SequenceUtils.nucleotides[r],mm[g][r][0],mm[g][r][1],mm[g][r][2]);
			dedupmm.close();
		}
		
		if (plot) {
			
			try {
				Logger.getGlobal().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".dedup.count.R");
				r.set("file",prefix+".dedup.count.tsv");
				r.set("out",prefix+".dedup.count.png");
				r.addSource(DedupUMI.class.getResourceAsStream("/resources/R/plotdedupcount.R"));
				r.run(false);

				r = new RRunner(prefix+".dedup.mm.R");
				r.set("file",prefix+".dedup.mm.tsv");
				r.set("out",prefix+".dedup.mm.png");
				r.addSource(DedupUMI.class.getResourceAsStream("/resources/R/plotdedupmm.R"));
				r.run(false);
				
				
				PlotReport pr = new PlotReport("Deduplication", StringUtils.toJavaIdentifier(inp+"_dedup"),
						FileUtils.getNameWithoutExtension(inp), "Distribution of deduplication events", prefix+".dedup.count.png", null, prefix+".dedup.count.R", prefix+".dedup.count.tsv");
				PlotReport pr2 = new PlotReport("Deduplication Mismatchcorrection", StringUtils.toJavaIdentifier(inp+"_dedup_Mismatches"),
						FileUtils.getNameWithoutExtension(inp), "Distribution of corrected mismatches", prefix+".dedup.mm.png", null, prefix+".dedup.mm.R", prefix+".dedup.mm.tsv");
				FileUtils.writeAllText(DynamicObject.from("plots",new Object[] {pr,pr2}).toJson(), new File(prefix+".dedup.report.json"));
				
			} catch (Throwable e) {
				Logger.getGlobal().log(Level.SEVERE, "Could not plot!", e);
			}
			
		}
		
	}
	
	private static ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> dedupUmi(
			ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read, IntArrayList dc,
			int[][][] mm, boolean nocollapse, boolean mismatchCorrection, int minReadCount) {
		
		return dedup(read, (data,c)->{
			HashMap<DnaSequence,ReadMerger> regions = new HashMap<>();
			for (int d=0; d<data.getDistinctSequences(); d++) {
//				if (data.getMultiplicity(d)!=1) throw new RuntimeException("Internal error: No multimappers allowed here!");
				
				for (DnaSequence bc : data.getBarcodes(d, c)) {
					ReadMerger p = regions.computeIfAbsent(bc, x->new ReadMerger());
					p.add(data,d, mm, 1);
				}
			}
			return regions;
		}, dc, mm, nocollapse, mismatchCorrection, minReadCount);
	}

	private static ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> dedupNoUmi(
			ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, IntArrayList dc,
			int[][][] mm, boolean nocollapse, boolean mismatchCorrection, int minReadCount) {
		
		return dedup(read, (data,c)->{
			HashMap<DnaSequence,ReadMerger> regions = new HashMap<>();
			DnaSequence bc = new DnaSequence();
			for (int d=0; d<data.getDistinctSequences(); d++) {
				int co = data.getCountFloor(d, c, ReadCountMode.Weight);
				if (co>0) {
					ReadMerger p = regions.computeIfAbsent(bc, x->new ReadMerger());
					p.multi=1;
					p.add(data,d, mm,co);
				}
			}
			return regions;
		}, dc, mm, nocollapse, mismatchCorrection, minReadCount);
	}
	
	private static <T extends AlignedReadsData> ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> dedup(
			ImmutableReferenceGenomicRegion<T> r, BiFunction<T, Integer, HashMap<DnaSequence,ReadMerger>> readToMergers, IntArrayList dc,
			int[][][] mm, boolean nocollapse, boolean mismatchCorrection, int minReadCount) {
		
		if (nocollapse) {
			AlignedReadsDataFactory fac = new AlignedReadsDataFactory(r.getData().getNumConditions());
			fac.start();
			fac.add(r.getData());
			return new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), fac.create());
		}
		
		ArrayList<DefaultAlignedReadsData> pre = new ArrayList<>();
		int[] expectedTotals = new int[r.getData().getNumConditions()];
		
//		if (r.toLocationString().equals("1+:2377603-2377753"))
//			System.out.println();
		
		for (int c=0; c<r.getData().getNumConditions(); c++) {
			HashMap<DnaSequence,ReadMerger> regions = readToMergers.apply(r.getData(), c);
			if (mismatchCorrection)
				correctMismatches(regions);
			
//			if (checkBcMM(regions)) {
//			}
			
			expectedTotals[c]=regions.size();
			
			for (DnaSequence pp : regions.keySet()) {
				
				ReadMerger p = regions.get(pp);
				
				dc.increment(p.count);

				if (p.count>=minReadCount) {
					AlignedReadsDataFactory fac = new AlignedReadsDataFactory(r.getData().getNumConditions());
					fac.start();
					fac.newDistinctSequence();
					fac.setMultiplicity(p.multi);
					fac.setCount(c, 1);
	//				if (p.geom!=null)
	//					fac.setGeometry(p.geom[0], p.geom[1], p.geom[2]);
					
					p.addVariationsAndGeom(fac,mm);
					
					
					
					DefaultAlignedReadsData ard = fac.create();
					pre.add(ard);
					if (pre.size()>10) {
						fac = new AlignedReadsDataFactory(r.getData().getNumConditions());
						fac.start();
						for (DefaultAlignedReadsData lard : pre) 
							fac.add(lard);
						fac.makeDistinct();
						pre.clear();
						pre.add(fac.create());
					}
				}
				else {
					expectedTotals[c]--;
				}
			}
		}
		
		if (pre.size()==0) {
			for (int i=0; i<expectedTotals.length; i++)
				if (expectedTotals[i]!=0)
					throw new RuntimeException("Internal error: Counts do not match!");
			return null;
		}
		
		ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> re;
		if (pre.size()==1)
			re = new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), pre.get(0));
		else {
			AlignedReadsDataFactory fac = new AlignedReadsDataFactory(r.getData().getNumConditions());
			fac.start();
			for (DefaultAlignedReadsData ard : pre) 
				fac.add(ard);
			fac.makeDistinct();
			re = new ImmutableReferenceGenomicRegion<>(r.getReference(),r.getRegion(),fac.create());
		}
		
		for (int i=0; i<expectedTotals.length; i++)
			if (expectedTotals[i]!=re.getData().getTotalCountForConditionInt(i, ReadCountMode.All))
				throw new RuntimeException("Internal error: Counts do not match!");
		
		return re;
	}
	
	
	private static boolean checkBcMM(HashMap<DnaSequence, ReadMerger> regions) {
		boolean re = false;
		DnaSequence[] a = regions.keySet().toArray(new DnaSequence[0]);
		for (int i=0; i<a.length; i++) {
			for (int j=i+1; j<a.length; j++) {
				if (StringUtils.hamming(a[i], a[j])<=1 && (regions.get(a[j]).count>=2*regions.get(a[i]).count-1||regions.get(a[i]).count>=2*regions.get(a[j]).count-1)) {
					System.out.println(a[i]+" "+regions.get(a[i]).count+"  "+a[j]+" "+regions.get(a[j]).count);
					re = true;
				}
			}
		}
		return re;
	}

	private static void correctMismatches(HashMap<DnaSequence, ReadMerger> regions) {
		
		// Do it like UMItools. We do it the other way around
		// I.e. each umi a is connected to an umi b with a directed edge a->b iff hamming distance is 1, and b has at least twice as many reads as a 
		// (-1, i.e. if both have 1, they are connected to each other) 
		// Then, each node with an outgoing edge must be added to its successor (the one with more reads, if there is more than one)
		
		// This does not contain the singletons!
		SimpleDirectedGraph<DnaSequence> gb = new MismatchGraphBuilder<DnaSequence>(regions.keySet().toArray(new DnaSequence[0]))
				.setInclude((a,b)->{
					return regions.get(b).count>=2*regions.get(a).count-1; 
				})
				.build();
			
		gb.sortAdjacencyLists((a,b)->Integer.compare(regions.get(b).count,regions.get(a).count));
		
//		for (DnaSequence s : gb.getTopologicalOrder()) {
//			if (gb.getTargets(s)!=null) {// this is sorted such that the first one is the one with the most reads
//				regions.get(gb.getTargets(s).node).add(regions.get(s));
//				regions.remove(s); // as only these barcodes are removed, singletons are retained!
//			}
//		}
		
		// Algorithm: start with the node with the fewest reads, and add it to the first successor (that is still there)
		DnaSequence[] nodes = regions.keySet().toArray(new DnaSequence[0]);
		Arrays.sort(nodes,(a,b)->Integer.compare(regions.get(a).count, regions.get(b).count));
		for (DnaSequence s : nodes) {
			if (!regions.containsKey(s)) 
				continue;
			
			// this is sorted such that the first one is the one with the most reads
			for (AdjacencyNode<DnaSequence> t = gb.getTargets(s); t!=null; t = t.next) {
				if (!regions.containsKey(t.node)) 
					continue;
				regions.get(t.node).add(regions.get(s));
				regions.remove(s); // as only these barcodes are removed, singletons are retained!
				break;
			}
		
		}

	}


	
	public static <T> void correctMismatchesSet(HashMap<String, HashSet<T>> regions) {
		
		// Do it like UMItools. We do it the other way around
		// I.e. each umi a is connected to an umi b with a directed edge a->b iff hamming distance is 1
		// (-1, i.e. if both have 1, they are connected to each other) 
		// Then, each node with an outgoing edge must be added to its successor (the one with more reads, if there is more than one)
		
		// This does not contain the singletons!
		SimpleDirectedGraph<String> gb = new MismatchGraphBuilder<String>(regions.keySet().toArray(new String[0]))
				.build();
			
		gb.sortAdjacencyLists((a,b)->Integer.compare(regions.get(b).size(),regions.get(a).size()));
		
		// Algorithm: start with the node with the fewest reads, and add it to the first successor (that is still there)
		String[] nodes = regions.keySet().toArray(new String[0]);
		Arrays.sort(nodes,(a,b)->Integer.compare(regions.get(a).size(), regions.get(b).size()));
		for (String s : nodes) {
			if (!regions.containsKey(s)) 
				continue;
			
			// this is sorted such that the first one is the one with the most reads
			for (AdjacencyNode<String> t = gb.getTargets(s); t!=null; t = t.next) {
				if (!regions.containsKey(t.node)) 
					continue;
				regions.get(t.node).addAll(regions.get(s));
				regions.remove(s); // as only these barcodes are removed, singletons are retained!
				break;
			}
		
		}

	}

	
	
	private static class ReadMerger extends HashMap<Integer,Counter<AlignedReadsMismatch>>{

		private ArrayList<VarsGeom> indelsAndGeom = new ArrayList<>();
		
		
		private int count = 0;
		private int multi = 99;
		
		@Override
		public String toString() {
			return ""+count;
		}
		
		public void addVariationsAndGeom(AlignedReadsDataFactory fac, int[][][] mm) {
			IntArrayList mms = new IntArrayList();
			ArrayList<AlignedReadsMismatch> toAdd = new ArrayList<AlignedReadsMismatch>(); // need this, because geom is determined afterwards, and first/second mismatches must be added!
			// actually, the size of indels should be the number of added reads or zero (if all reads had the same length), and the cumulated length deltas by insertion or deletion should be the same for all elements (very rarely: read 1 might have no indel, read 2 might have insertion of length l and a deletion of length l) 
			for (Integer pos : keySet()) {
				// majority vote
				Counter<AlignedReadsMismatch> ctr = get(pos);
				AlignedReadsMismatch maxVar = ctr.getMaxElement(0);
				int total = ctr.total()[0];
				if (maxVar!=null && ctr.get(maxVar)[0]>count-total) {
					toAdd.add(maxVar);
					getForVar(mm, maxVar)[0]++;
					getForVar(mm, maxVar)[1]+=ctr.get(maxVar, 0);
					mms.add(pos);
				}
			}
			
			
			// decide on indel pattern!
			VarsGeom max = isAllTheSame(indelsAndGeom);
			
			if (max==null) {
				HashMap<VarsGeom,MutableInteger> counter = new HashMap<VarsGeom, MutableInteger>();
				int maxval = 0;
				
				for (VarsGeom pair : indelsAndGeom) {
					AlignedReadsVariation[] pattern = pair.indels;
					// check if compatible with mms
					int numMMinPattern = EI.wrap(pattern).filter(v->v.isDeletion() && findPosInDeletion(mms,v.asDeletion())).countInt();
					if (numMMinPattern==0) {
						// we only get here if the pattern is compatible: count it for majority vote
						MutableInteger cc = counter.computeIfAbsent(pair, x->new MutableInteger());
						cc.N++;
						if (cc.N>maxval) {
							maxval=cc.N;
							max=pair;
						}
					}
				}
			}
			
			if (max==null) {
				// we only get here if there is no compatible pattern
				// there is nothing we can do at this point (other than removing mismatches, but hey...)
				// so we don't do anything and live with the fact that the read now has a different length
				max=indelsAndGeom.get(0);
				max.indels = new AlignedReadsVariation[0];
			}
			
			// add softclips from all indelsAndGeoms that are equal to max (i.e. have the same indel pattern and softclip lengths)
			VarsGeom umax = max;
			fac.addSoftclip(true, getSoftclip(EI.wrap(indelsAndGeom).filter(r->r.equals(umax)).map(r->r.softclips[0])),false);
			fac.addSoftclip(false, getSoftclip(EI.wrap(indelsAndGeom).filter(r->r.equals(umax)).map(r->r.softclips[1])),false);
			fac.addSoftclip(true, getSoftclip(EI.wrap(indelsAndGeom).filter(r->r.equals(umax)).map(r->r.softclips[2])),true);
			fac.addSoftclip(false, getSoftclip(EI.wrap(indelsAndGeom).filter(r->r.equals(umax)).map(r->r.softclips[3])),true);
			
//			addSoft(fac,true,s5p[0],false);
//			addSoft(fac,true,s5p[1],true);
//			addSoft(fac,false,s3p[0],false);
//			addSoft(fac,false,s3p[1],true);

			
			int[] geom = max.geom;
			// now finally add the mismatches
			for (AlignedReadsMismatch mmv : toAdd) {
				fac.addVariation(mmv);
				if (geom!=null && mmv.getPosition()>=geom[0] && mmv.getPosition()<geom[0]+geom[1]) 
					fac.addVariationToOtherRead(mmv);
			}
			
			
			for (AlignedReadsVariation v : max.indels)
				fac.addVariation(v);
			if (geom!=null)
				fac.setGeometry(geom[0], geom[1], geom[2]);
				
		}
		
		private VarsGeom isAllTheSame(
				ArrayList<VarsGeom> l) {
			for (int i=1; i<l.size(); i++)
				if (!l.get(0).equals(l.get(i)))
					return null;
			return l.get(0);
		}

		private boolean findPosInDeletion(IntArrayList pos, AlignedReadsDeletion dele) {
			// pos might not be sorted, so binary search is out of question!
			return pos.iterator().filterInt(p->dele.contains(p)).sum()>0;
		}

		private CharSequence getSoftclip(ExtendedIterator<AlignedReadsSoftclip> softclips) {
			ArrayList<CharSequence> l = softclips.map(s->s.getReadSequence()).list();
//		private void addSoft(AlignedReadsDataFactory fac, boolean s5p, ArrayList<CharSequence> l, boolean secondRead) {
			if (l.size()==0) 
				return "";
			if (l.size()==1 || EI.wrap(l).skip(1).filter(c->l.get(0).equals(c)).count()==l.size()-1) {
//				fac.addSoftclip(s5p, l.get(0), secondRead);
				return l.get(0);
			}
			else {
				char[] m = new char[EI.wrap(l).mapToInt(c->c.length()).unique(true).getUniqueResult(true, true)];
				int[] co = new int[5];
				for (int p=0; p<m.length; p++) {
					Arrays.fill(co, 0);
					for (CharSequence c : l) 
						if (p<c.length()) {
							int b = SequenceUtils.inv_nucleotides[c.charAt(p)];
							if (b>=0 && b<4)
								co[b]++;
							else
								co[4]++;
						}
					int cop = ArrayUtils.argmax(co);
					m[p] = SequenceUtils.nucleotides[cop];
				}
				return String.valueOf(m);
//				fac.addSoftclip(s5p, String.valueOf(m, 0, len), secondRead);
			}
		}

		private int[] getForVar(int[][][] mm, AlignedReadsVariation var) {
			return mm[SequenceUtils.inv_nucleotides[var.getReferenceSequence().charAt(0)]][SequenceUtils.inv_nucleotides[var.getReadSequence().charAt(0)]];
		}
		
		public void add(ReadMerger other) {
			count+=other.count;
			multi = Math.min(multi, other.multi);
			indelsAndGeom.addAll(other.indelsAndGeom);
			for (Integer pos : other.keySet())
				computeIfAbsent(pos, x->new Counter<>()).integrate(other.get(pos));
		}

		public void add(AlignedReadsData r, int d, int[][][] mm, int num) {
			count+=num;
			multi = Math.min(multi, r.getMultiplicity(d));
			int[] geom = r.hasGeometry()?new int[] {r.getGeometryBeforeOverlap(d),r.getGeometryOverlap(d),r.getGeometryAfterOverlap(d)}:null;
			
			ArrayList<AlignedReadsVariation> indel = new ArrayList<AlignedReadsVariation>(); 
			ArrayList<AlignedReadsSoftclip> soft = new ArrayList<AlignedReadsSoftclip>(); 
			int vc = r.getVariationCount(d);
			for (int v=0; v<vc; v++) {
				AlignedReadsVariation vari = r.getVariation(d, v);
				
				boolean consistent = true;
				boolean overlap = r.hasGeometry() && r.isPositionInOverlap(d, vari.getPosition());
				
				// remove inconsistent variations (first vs second read)
				if (overlap) {
					int offset = vari.isFromSecondRead()?-1:1;
					if (v+offset>=vc || v+offset<0 || !vari.isConsistentInOtherRead(r.getVariation(d, v+offset)))
						consistent=false;
				}
				if (vari.isMismatch()) {
					int pos = vari.getPosition();
					if (!overlap || (!vari.isFromSecondRead() && consistent)) // only add first read mismatches
						computeIfAbsent(pos, x->new Counter<>()).add(vari.asMismatch());
					getForVar(mm, vari)[2]++;
				} else if (vari.isSoftclip()) {
					soft.add(vari.asSoftclip());
				} else if (vari.isDeletion() || vari.isInsertion())
					indel.add(vari);
			}
			VarsGeom vg = new VarsGeom(indel.toArray(new AlignedReadsVariation[0]),soft.toArray(new AlignedReadsSoftclip[0]),geom);
			for (int i=0; i<num; i++)
				indelsAndGeom.add(vg);
		}
		
	}
	
	/*
	 * Does only compare softclip *lengths*
	 */
	private static AlignedReadsSoftclip[] emptySoftclips = {
			new AlignedReadsSoftclip(true, "", false),	
			new AlignedReadsSoftclip(false, "", false),	
			new AlignedReadsSoftclip(true, "", true),	
			new AlignedReadsSoftclip(false, "", true)
	};
	private static class VarsGeom {
		AlignedReadsVariation[] indels;
		AlignedReadsSoftclip[] softclips;  // 5p,3p,5pr,3pr
		int[] geom;
		public VarsGeom(AlignedReadsVariation[] indels, AlignedReadsSoftclip[] softclips, int[] geom) {
			super();
			this.indels = indels;
			this.softclips = new AlignedReadsSoftclip[4];
			for (AlignedReadsSoftclip s : softclips)
				this.softclips[ (s.isFromSecondRead()?2:0) + (s.getPosition()==0?0:1) ] = s;
			for (int i=0; i<softclips.length; i++)
				if (softclips[i]==null)
					softclips[i]=emptySoftclips[i];
			this.geom = geom;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(geom);
			result = prime * result + Arrays.hashCode(indels);
			result = prime * result + lenHashCode(softclips);
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			VarsGeom other = (VarsGeom) obj;
			if (!Arrays.equals(geom, other.geom))
				return false;
			if (!Arrays.equals(indels, other.indels))
				return false;
			if (!lenEquals(softclips, other.softclips))
				return false;
			return true;
		}
		private static int lenHashCode(AlignedReadsSoftclip[] a) {
	        if (a == null)
	            return 0;

	        int result = 1;

	        for (AlignedReadsSoftclip element : a)
	            result = 31 * result + (element == null ? 0 : element.getReadSequence().length());

	        return result;
	    }
		private static boolean lenEquals(AlignedReadsSoftclip[] a,AlignedReadsSoftclip[] a2) {
			if (a==a2)
	            return true;
	        if (a==null || a2==null)
	            return false;

	        int length = a.length;
	        if (a2.length != length)
	            return false;

	        for (int i=0; i<length; i++) {
	        	AlignedReadsSoftclip o1 = a[i];
	        	AlignedReadsSoftclip o2 = a2[i];
	            if (!(o1==null ? o2==null : (o2!=null && o1.getReadSequence().length()==o2.getReadSequence().length())))
	                return false;
	        }

	        return true;
	    }
		@Override
		public String toString() {
			return "VarsGeom [indels=" + Arrays.toString(indels) + ", softclips=" + Arrays.toString(softclips)
					+ ", geom=" + Arrays.toString(geom) + "]";
		}
		
		
	}
}
