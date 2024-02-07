package executables;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.TreeSet;

import cern.colt.bitvector.BitVector;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.SpliceGraph;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.io.text.HeaderLine;
import gedi.util.math.stat.counting.Counter;
import gedi.util.math.stat.descriptive.WeightedMeanVarianceOnline;
import gedi.util.math.stat.kernel.EpanechnikovKernel;
import gedi.util.math.stat.kernel.PreparedIntKernel;
import gedi.util.mutable.MutableInteger;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;
import gedi.util.sequence.DnaSequence;
import jdistlib.Normal;


public class Read10XMem {
	
	




	public static void main(String[] args) throws IOException {
		
		
		Read10XParameterSet params = new Read10XParameterSet();
		GediProgram pipeline = GediProgram.create("Read10X",
				new Read10XInfer3pProgram(params),
				new Read10XClusterProgram(params)
				);
		GediProgram.run(pipeline, null, null, new CommandLineHandler("Read10X","Read10X reads 10x bam files and creates a demultiplexed and de-barcoded cit file.",args));

	}
	
	public enum PeakCategory {
		Exonic,Spliced,ExonIntronic,Intronic,Flank,Intergenic
	}

	
	public static ArrayList<GenomicRegion> getFirstPolyAStretch(String seq) {
		
		ArrayList<GenomicRegion> re = new ArrayList<>();
		
		// A/G stretch in 60-120? (len 10, >50%A)
		int[] cumA = new int[seq.length()];
		int[] cumG = new int[seq.length()];
		for (int i=0; i<seq.length(); i++) {
			if (seq.charAt(i)=='A')
				cumA[i]++;
			else if (seq.charAt(i)=='G')
				cumG[i]++;
		}
		ArrayUtils.cumSumInPlace(cumA, 1);
		ArrayUtils.cumSumInPlace(cumG, 1);
		for (int i=10; i<cumA.length; i++) {
			int na=cumA[i]-cumA[i-10];
			int ng=cumG[i]-cumG[i-10];
			if (na+ng>=9 && na>ng) {
				// extend
				int s;
				for (s=i+1; s<cumA.length; s++) {
					na=cumA[s]-cumA[i-10];
					ng=cumG[s]-cumG[i-10];
					if (na+ng<s-i+9 || na<=ng) 
						break;
				}
				re.add(new ArrayGenomicRegion(i-10,s));
			}
		}
		
		return re;
	}
	
	
	
	public static class Read10XClusterProgram extends GediProgram {



		public Read10XClusterProgram(Read10XParameterSet params) {
			addInput(params.prefix);
			addInput(params.reads);
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.maxdist);
			addInput(params.whitelist);
			addInput(params.test);
			addInput(params.infer3pparam);

			addOutput(params.clustercit);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			GenomicRegionStorage<BarcodedAlignedReadsData> reads = getParameter(1);
			Genomic genomic = getParameter(2);
			int nthreads = getIntParameter(3);
			int maxdist = getIntParameter(4);
			String whitelistFile = getParameter(5);
			boolean test = getBooleanParameter(6);
			File pparam = getParameter(7);

			
			double[] meansd = {Double.NaN, Double.NaN};
			HeaderLine h = new HeaderLine();
			EI.lines(pparam).header(h).split("\t").forEachRemaining(a->{
				if (a[0].equals("Mean"))
					meansd[0] = Double.parseDouble(a[1]);
				else if (a[0].equals("Sd"))
					meansd[1] = Double.parseDouble(a[1]);
			});
			
			context.logf("Mean=%.1f, Sd=%.1f", meansd[0],meansd[1]);
			
			PreparedIntKernel kernel = new EpanechnikovKernel(meansd[1]/4).prepare(true);

			
			genomic.getTranscriptTable("source");
			
			
			context.getLog().info("Identifying clusters...");
			
			String[] conditions = reads.getMetaDataConditions();
			HashMap<String, Integer> icond = EI.wrap(conditions).indexPosition();
			
			HashMap<String,Integer>[] whitelists = new HashMap[conditions.length];
			for (int i=0; i<whitelists.length; i++)
				whitelists[i] = new HashMap<>();
			
			ArrayList<DynamicObject> meta = new ArrayList<>();
			h = new HeaderLine();
			int ncond = 0;
			for (String[] a : EI.lines(whitelistFile).header(h).split('\t').loop()) {
				Integer ind = icond.get(a[h.get("Condition")]);
				if (ind==null) throw new RuntimeException("Condition "+a[h.get("Condition")]+" as given in the whitelist file is unknown!");
				whitelists[ind].put(a[h.get("Barcode")], ncond++);
				meta.add(DynamicObject.from("name", a[h.get("Condition")]+"."+a[h.get("Barcode")]));
			}
			
			TreeSet<ReferenceSequence> refs = new TreeSet<ReferenceSequence>(reads.getReferenceSequences());
			refs.removeIf(r->!genomic.getSequenceNames().contains(r.getName()));
			
			CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(getOutputFile(0).getPath(),DefaultAlignedReadsData.class);
			out.setMetaData(DynamicObject.from("conditions",meta));
			
//			MutableInteger id = new MutableInteger();
			int block = 1;
			IterateIntoSink<ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>>> sink = new IterateIntoSink<>(
					ei->out.fill(ei.parallelized(nthreads, block, xei->
									xei.unfold(x->{
										ImmutableReferenceGenomicRegion<Void>[] peaks = findPeaks(x, conditions, meansd[0], meansd[1], kernel, genomic);
										if (peaks.length==0) return null;
										return demultiplex(x, peaks, meansd[0],meansd[1],conditions, genomic, whitelists);
//										return EI.<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>empty();
									}
											))
										.removeNulls()
										.progress(context.getProgress(), -1, x->"Detecting clusters "+x.toLocationString())
										,context.getProgress()),
//					ei->out.fill(ei.map(r->new ImmutableReferenceGenomicRegion<NameAnnotation>(r.getReference(), r.getRegion(),new NameAnnotation(id.N+++"")))),
					nthreads*block
					);
			
			
			IntervalTree<GenomicRegion,ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> openClusters = null;
			TreeMap<Integer,GenomicRegion> openClustersByStop = null;
			
			MemoryIntervalTreeStorage<Transcript> transcripts = genomic.getTranscripts();
			TrimmedGenomicRegion reuse = new TrimmedGenomicRegion();
			
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : (test?reads.ei("JN555585+"):reads.ei()) // reads.ei("10+:987573-1020065").sort() //10+:1043717-1044348
										.progress(context.getProgress(), (int)reads.size(), r->"Collecting reads "+r.toLocationString())
										.filter(r->r.getData().getMultiplicity(0)<=1)
										.checkOrder((a,b)->a.compareTo(b))
										.loop()) {
				
				ArrayList<ImmutableReferenceGenomicRegion<Transcript>> ctrans = transcripts.ei(read).chain(transcripts.ei(read).map(t->t.toMutable().toOppositeStrand().toImmutable())).filter(t->isCompatible(t.getRegion(),reuse.set(read.getRegion()))).list();
				
				if (read.getRegion().getNumParts()>1 && ctrans.isEmpty())
					continue;
				
				ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>> re;
				if (openClusters==null) {
					// first cluster
					openClusters = new IntervalTree<GenomicRegion, ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>>(read.getReference());
					openClusters.put(read.getRegion(),re = new ArrayList<>());
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				}
				else if (!openClusters.getReference().equals(read.getReference())) {
					// first cluster on this chromosome
					
					for (GenomicRegion clreg : openClusters.keySet()) 
						sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, openClusters.get(clreg)));

					openClusters = new IntervalTree<GenomicRegion, ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>>(read.getReference());
					openClusters.put(read.getRegion(),re = new ArrayList<>());
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				} else {
					ArrayList<GenomicRegion> overlaps = openClusters
						.keys(read.getRegion().getStart(), read.getRegion().getStop())
						.filter(reg->intersects(reuse.set(reg),read.getRegion(),ctrans,maxdist)).list();
					
					if (overlaps.size()>0) {
						// merge clusters 
						re = openClusters.get(overlaps.get(0));
						GenomicRegion r = overlaps.get(0);
						for (int i=1; i<overlaps.size(); i++) {
							re.addAll(openClusters.get(overlaps.get(i)));
							r = r.union(overlaps.get(i));
						}
						r = r.union(read.getRegion());
						
						for (GenomicRegion reg : overlaps) {
							openClusters.remove(reg);
							openClustersByStop.remove(reg.getStop());
						}
						
						openClusters.put(r, re);
						openClustersByStop.put(r.getStop(), r);

					} else {
						openClusters.put(read.getRegion(), re = new ArrayList<>());
						openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
					}
					
					// prune interval tree (everything that ends left of the current start
					NavigableMap<Integer, GenomicRegion> head = openClustersByStop.headMap(read.getRegion().getStart(), false);
					for (GenomicRegion clreg : head.values()) 
						sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, openClusters.get(clreg)));

					openClusters.keySet().removeAll(head.values());
					head.clear();
				}
				re.add(read);
				
			}
			
			for (GenomicRegion clreg : openClusters.keySet()) 
				sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, openClusters.get(clreg)));

			
			sink.finish();
			

			return null;
		}
		
		
		

		@SuppressWarnings("unchecked")
		private static ImmutableReferenceGenomicRegion<Void>[] findPeaks(
				ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> or, 
				String[] conditions, double mean, double sd, PreparedIntKernel kernel, Genomic g) {
			
			System.out.println("Cluster "+or.toLocationString()+" "+or.getData().size());
			
			MutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> r = or.toMutable();
			
			BarcodeCounter[] bco = new BarcodeCounter[conditions.length];
			for (int c=0; c<conditions.length; c++) 
				bco[c] = new BarcodeCounter();

			// determine weight for each barcode
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : r.getData()) {
				for (int c=0; c<conditions.length; c++)
					for (int d=0; d<read.getData().getDistinctSequences(); d++)
						bco[c].addBarcodes(read.getData().getBarcodes(d, c));
			}
			System.out.println("Determined barcodes: "+EI.wrap(bco).mapToInt(bc->bc.counter.size()).sum());
			
			HashMap<DnaSequence,MutableInteger>[] w = new HashMap[bco.length];
			for (int c=0; c<bco.length; c++)
				w[c] = bco[c].get();
			// If a umi occurs n times, its weight is n (such that it can be counted as 1/n)
			
			// map reads to their weights
			ArrayList<ImmutableReferenceGenomicRegion<double[]>> umiw = readsToWeights(r.getData(),w);
			
			// kernel density estimation, alters r!
			double[][] smoothed = smooth(g,r,umiw,sd);

			System.out.println("smoothed");
//			try {
//				System.out.println(r.toLocationString());
//				EI.wrap(smoothed[0]).print("smoothed");
//				EI.wrap(smoothed[1]).print("smoothed2");
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
			
//			kernel.processInPlace(smoothed[0], 0, smoothed[0].length);
//			
//			System.out.println("smoothed 2");
//			try {
//				EI.wrap(smoothed[0]).print("ssmoothed");
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
			
			
			// identify peaks
			int[] peaks = detectPeaks(g,smoothed,r, mean, sd);
			ImmutableReferenceGenomicRegion<Void>[] re = EI.wrap(peaks)//.filter(p->or.getRegion().contains(r.map(new ArrayGenomicRegion(p,p+1))))
						.map(p->new ImmutableReferenceGenomicRegion<Void>(r.getReference(),r.mapMaybeOutSide(new ArrayGenomicRegion(p,p+1))))
						.toArray(new ImmutableReferenceGenomicRegion[0]);
			
			System.out.println("Peaks: "+peaks.length);
//			try {
//				EI.wrap(peaks).print("peaks");
//			} catch (IOException e) {
//				e.printStackTrace();
//			}

			return re;
		}

		private static int[] detectPeaks(Genomic g, 
				double[][] smoothed12,
				MutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> r, double mean, double bw) {
			
			double[] smoothed = smoothed12[0];
			double[] smoothed2 = smoothed12[1];
			
			int maxo = (int) Normal.quantile(0.99, 0, bw, true, false);
			
			SpliceGraph introns = new SpliceGraph(0,r.getRegion().getTotalLength());
			g.getTranscripts().ei(r).chain(g.getTranscripts().ei(r).map(t->t.toMutable().toOppositeStrand().toImmutable()))
			.forEachRemaining(t->introns.addIntrons(r.induce(t.getRegion())));
		
			IntArrayList re = new IntArrayList();
			
			for (;;) {
				
				int argmax = ArrayUtils.argmax(smoothed);
				if (smoothed[argmax]==0) {
					re.sort();
					re.unique();
					return re.toIntArray();
				}
//				System.out.println(argmax+" "+smoothed[argmax]);
				
				LinkedList<ProfileWalker> walkers = new LinkedList<>();
				ArrayList<ProfileWalker> add = new ArrayList<>();
				walkers.add(new ProfileWalker(argmax,-1,smoothed,null));
				walkers.add(new ProfileWalker(argmax,1,smoothed,smoothed2));
				smoothed[argmax] = 0;
				
				while (!walkers.isEmpty()) {
					
					Iterator<ProfileWalker> it = walkers.iterator();
					while (it.hasNext()) {
						ProfileWalker w = it.next();
						if (!w.moveAndCheckIntrons(introns, add,maxo)) {
							it.remove();
							if (w.bestDiffPos>=0 && w.hasPathRead()) {
								re.add(w.bestDiffPos);
//								System.out.println(" "+w.bestDiffPos);
							}
							
						}
					}
					walkers.addAll(add);
					add.clear();
				}
				
//				try {
//					EI.wrap(smoothed).print("smoothed");
//				} catch (IOException e) {
//					e.printStackTrace();
//				}
				
//				IntArrayList ppu = new IntArrayList();
//				DoubleArrayList oldvu = new DoubleArrayList();
//				ppu.add(argmax);
//				oldvu.add(smoothed[argmax]);
//				IntArrayList ppd = new IntArrayList();
//				IntArrayList stopcand = new IntArrayList();
//				IntArrayList ppdex = new IntArrayList();
//				DoubleArrayList oldvd = new DoubleArrayList();
//				DoubleArrayList oldvdex = new DoubleArrayList();
//				ppd.add(argmax);
//				stopcand.add(-1);
//				oldvd.add(smoothed[argmax]);
//				
//				HashSet<Intron> usedIntron = new HashSet<>();
//				
//				int o = 1;
//				// upstream
//				while (!ppu.isEmpty() || !ppd.isEmpty()) {
//					int uo = o;
//					for (int ind=0; ind<ppu.size(); ind++) {
//						int ee = ppu.getInt(ind)-o;
//						double oldvv = oldvu.getDouble(ind);
//						introns.getIntrons().forEachIntervalIntersecting(ppu.getInt(ind)-o, ppu.getInt(ind)-o, intron->{
//							if (intron.getStop()==ee && usedIntron.add(intron)) {
//								ppu.add(intron.getStart()-1+uo);
//								oldvu.add(oldvv);
//								ppdex.add(intron.getStart()-uo);
//								oldvdex.add(oldvv);
//							}
//						});
//					}
//					for (int ind=0; ind<ppd.size(); ind++) {
//						int ee = ppd.getInt(ind)+o;
//						double oldvv = oldvd.getDouble(ind);
//						introns.getIntrons().forEachIntervalIntersecting(ppd.getInt(ind)+o, ppd.getInt(ind)+o, intron->{
//							if (intron.getStart()==ee && usedIntron.add(intron)) {
//								ppd.add(intron.getEnd()-uo);
//								stopcand.add(-1);
//								oldvd.add(oldvv);
//								ppu.add(intron.getEnd()-1+uo);
//								oldvu.add(oldvv);
//							}
//						});
//					}
//					for (int ind=0; ind<ppu.size(); ind++)
//						if (ppu.getInt(ind)-o>=0 && (smoothed[ppu.getInt(ind)-o]<oldvu.getDouble(ind) || o<maxo)) {
//							oldvu.set(ind, smoothed[ppu.getInt(ind)-o]);
//							smoothed[ppu.getInt(ind)-o]=0;
//						}
//						else {
//							oldvu.removeEntry(ind);
//							ppu.removeEntry(ind--);
//						}
//					for  (int ind=0; ind<ppd.size(); ind++) {
//						// check last read
//						if (ppd.getInt(ind)+o<smoothed2.length && smoothed2[ppd.getInt(ind)+o]<smoothed2[ppd.getInt(ind)+o-1]) {
//							int sc = stopcand.getInt(ind);
//							if (sc==-1) {
//								stopcand.set(ind,ppd.getInt(ind)+o);
//							}
//							else {
//								double ldiff = smoothed2[sc-1]-smoothed2[sc];
//								double hdiff = smoothed2[ppd.getInt(ind)+o-1]-smoothed2[ppd.getInt(ind)+o];
//								
//								if (hdiff>=0.5*ldiff)
//									stopcand.set(ind,ppd.getInt(ind)+o);
//							}
//							
//						}
//
//						if (ppd.getInt(ind)+o<smoothed.length && (smoothed[ppd.getInt(ind)+o]<oldvd.getDouble(ind) || o<maxo)) {
//							oldvd.set(ind, smoothed[ppd.getInt(ind)+o]);
//							smoothed[ppd.getInt(ind)+o]=0;
//						}
//						else {
//							if (smoothed2[argmax]>0) {
//								int sc = stopcand.getInt(ind);
//								if (sc==-1) {
//									re.add(ppd.getInt(ind)+Math.max((int)mean, o));
//									System.out.println(" "+ind+" "+(ppd.getInt(ind)+Math.max((int)mean, o)));
//								} else {
//									re.add(sc+Math.max((int)mean-o, 0));
//									System.out.println(" "+ind+" "+(sc+Math.max((int)mean-o, 0)));
//								}
//							}
//
//							oldvd.removeEntry(ind);
//							ppd.removeEntry(ind--);
//						}
//					}
//					for (int ind=0; ind<ppdex.size(); ind++)
//						if (ppdex.getInt(ind)+o<smoothed.length && (smoothed[ppdex.getInt(ind)+o]<oldvdex.getDouble(ind) || o<maxo)) {
//							oldvdex.set(ind, smoothed[ppdex.getInt(ind)+o]);
//							smoothed[ppdex.getInt(ind)+o]=0;
//						}
//						else {
//							oldvdex.removeEntry(ind);
//							ppdex.removeEntry(ind--);
//						}
//					o++;
//				}
//				smoothed[argmax] = 0;
				
			}
			
		}

		// [density,density restricted to read cover]
		private static double[][] smooth(Genomic g, MutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> r, ArrayList<ImmutableReferenceGenomicRegion<double[]>> umis, double bw) {
			
			int maxo = (int) Normal.quantile(0.99, 0, bw, true, false);
			
			r.alterRegion(re->re.extendAll(maxo, maxo));
			SpliceGraph introns = new SpliceGraph(0,r.getRegion().getTotalLength());
			g.getTranscripts().ei(r).chain(g.getTranscripts().ei(r).map(t->t.toMutable().toOppositeStrand().toImmutable()))
				.forEachRemaining(t->introns.addIntrons(r.induce(t.getRegion())));
			
			double[] re = new double[r.getRegion().getTotalLength()];
			double[] re2 = new double[r.getRegion().getTotalLength()];
			
			double[] w = new double[maxo];
			for (int i=0; i<maxo; i++)
				w[i] = Normal.density(i, 0, bw, false);
								
			for (ImmutableReferenceGenomicRegion<double[]> read : umis) {
				
				double rw = ArrayUtils.sum(read.getData());
				int p = GenomicRegionPosition.FivePrime.position(read);
				int mp = r.induce(p);
				
				IntArrayList pp = new IntArrayList();
				pp.add(mp);
				
				// upstream and pos itself
				for (int o=0; o<maxo & mp-o>=0; o++) {
					int uo = o;
					for (int ind=0; ind<pp.size(); ind++) {
						int ee = pp.getInt(ind)-o;
						introns.getIntrons().forEachIntervalIntersecting(pp.getInt(ind)-o, pp.getInt(ind)-o, intron->{
							if (intron.getStop()==ee)
								pp.add(intron.getStart()-1+uo);
						});
					}
					
					for (int ind=0; ind<pp.size(); ind++)
						if (pp.getInt(ind)-o>=0 && pp.getInt(ind)-o<re.length) {
							re[pp.getInt(ind)-o]+=w[o]*rw/pp.size();
						}
				}
				// down the read
//				re2[r.induce(read.map(0))]+=w[0]*rw;
				for (int o=1; o<Math.min(maxo, read.getRegion().getTotalLength()); o++) {
					int rp = r.induce(read.map(o));
					re[rp]+=w[o]*rw;
					if (o>2 && o<read.getRegion().getTotalLength()-3)
						re2[rp]+=w[o]*rw;
				}
				
				pp.clear();
				pp.add(r.induce(GenomicRegionPosition.ThreePrime.position(read, 1)));
				// downstream of the read
				for (int o=0; o<maxo-read.getRegion().getTotalLength(); o++) {
					int uo = o;
					for (int ind=0; ind<pp.size(); ind++) {
						int ee = pp.getInt(ind)+o;
						introns.getIntrons().forEachIntervalIntersecting(pp.getInt(ind)+o, pp.getInt(ind)+o, intron->{
							if (intron.getStart()==ee)
								pp.add(intron.getEnd()-uo);
						});
					}
					
					for (int ind=0; ind<pp.size(); ind++)
						if (pp.getInt(ind)+o<re.length && pp.getInt(ind)+o>=0) {
							re[pp.getInt(ind)+o]+=w[o+read.getRegion().getTotalLength()]*rw/pp.size();
						}
				}
			}
			
			return new double[][] {re,re2};
		}

		private static ArrayList<ImmutableReferenceGenomicRegion<double[]>> readsToWeights(
				ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>> data,
				HashMap<DnaSequence, MutableInteger>[] w) {
			ArrayList<ImmutableReferenceGenomicRegion<double[]>> re = new ArrayList<>();
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : data) {
				double[] a = new double[read.getData().getNumConditions()];
				for (int c=0; c<a.length; c++)
					for (int d=0; d<read.getData().getDistinctSequences(); d++)
						for (DnaSequence bc : read.getData().getBarcodes(d, c)) {
							a[c]+=1.0/w[c].get(bc).N;
						}
				re.add(new ImmutableReferenceGenomicRegion<>(read.getReference(), read.getRegion(),a));
			}
			return re;
		}
		
		private static ExtendedIterator<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> demultiplex(
				ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> or, 
				ImmutableReferenceGenomicRegion<Void>[] peaks, double mean, double sd, String[] conditions, 
				Genomic g,
				HashMap<String,Integer>[] whitelists) {
			
			TrimmedGenomicRegion reuse = new TrimmedGenomicRegion();
			Normal dist = new Normal(mean, sd);
			int cellBarcodeLength = whitelists[0].keySet().iterator().next().length();
			int nind = EI.wrap(whitelists).unfold(m->EI.wrap(m.values())).mapToInt(x->x).max()+1;
			
			// collect all reads per barcode (full reg, raw probs to peaks, em probs to peaks
			HashMap<MutablePair<DnaSequence,Integer>,MutableTriple<GenomicRegion,double[],Read10XMerger>> regions = new HashMap<>();
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : or.getData()) {
				for (int d=0; d<r.getData().getDistinctSequences(); d++) {
					for (int c=0; c<r.getData().getNumConditions(); c++) {
						for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
							String cellbc = bc.toString().substring(0,cellBarcodeLength);
							Integer cell = whitelists[c].get(cellbc); 
							if (cell!=null) {
								MutablePair<DnaSequence, Integer> key = new MutablePair<DnaSequence,Integer>(bc,c);
								MutableTriple<GenomicRegion,double[],Read10XMerger> p = regions.computeIfAbsent(key, x->new MutableTriple<>(new ArrayGenomicRegion(),new double[peaks.length],null));
								p.Item1 = p.Item1.union(r.getRegion());
							}
						}
					}
				}
			}
			
//			System.out.println("Collected reads: "+regions.size());
//			int n=0;
			
			Read10XEmObject em = new Read10XEmObject(peaks);
			Iterator<MutableTriple<GenomicRegion, double[], Read10XMerger>> it = regions.values().iterator();
			while (it.hasNext()) {
				MutableTriple<GenomicRegion, double[], Read10XMerger> rdp = it.next();
				ArrayList<GenomicRegion> compt = g.getTranscripts().ei(or.getReference(),rdp.Item1).chain(g.getTranscripts().ei(or.getReference().toOppositeStrand(),rdp.Item1))
						.filter(t->t.getRegion().contains(rdp.Item1))
						.map(t->t.getRegion())
						.list();
				compt.add(rdp.Item1.removeIntrons());
				Collections.sort(compt,(a,b)->{
					GenomicRegion ai = a.induce(rdp.Item1);
					GenomicRegion bi = b.induce(rdp.Item1);
					return Integer.compare(ai.getEnd()-ai.getStart(),bi.getEnd()-bi.getStart());
				});
				rdp.Item1 = compt.get(0).intersect(new ArrayGenomicRegion(rdp.Item1.getStart(),rdp.Item1.getEnd()));
				if (!em.add(g,dist,rdp,reuse))
					it.remove();
//				n++;
//				if (n%10000==0) System.out.println(n);
			}
			if (regions.isEmpty()) return EI.empty();
			
//			System.out.println("Start em");
			em.iterate(regions.values(),10);
			while (em.checkUnnecessary(regions.values(),0.1))
				em.iterate(regions.values(),10);
			
//			System.out.println("Finished em");
			
			HashMap<MutablePair<DnaSequence,Integer>,ImmutableReferenceGenomicRegion<Void>> barcodeToPeak = em.getMap(regions);
			
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : or.getData()) {
				for (int d=0; d<r.getData().getDistinctSequences(); d++) {
					for (int c=0; c<r.getData().getNumConditions(); c++) {
						for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
							String cellbc = bc.toString().substring(0,cellBarcodeLength);
							Integer cell = whitelists[c].get(cellbc); 
							if (cell!=null) {
								MutablePair<DnaSequence, Integer> key = new MutablePair<DnaSequence,Integer>(bc,c);
								MutableTriple<GenomicRegion,double[],Read10XMerger> p = regions.get(key);
								if (p!=null) { // because of removal
									if (p.Item3==null) {
										
										ImmutableReferenceGenomicRegion<Void> peak = barcodeToPeak.get(key);
										
										GenomicRegion reg = extend(g, p.Item1, peak, dist,reuse);
										ImmutableReferenceGenomicRegion<Void> rgr = new ImmutableReferenceGenomicRegion<>(or.getReference(), reg);
										if (rgr.getRegion().getTotalLength()>DefaultAlignedReadsData.MAX_POSITION) {
											reg = rgr.map(new ArrayGenomicRegion(rgr.getRegion().getTotalLength()-DefaultAlignedReadsData.MAX_POSITION,rgr.getRegion().getTotalLength()));
											rgr = new ImmutableReferenceGenomicRegion<>(or.getReference(), reg);
										}
										p.Item1 = reg;
										p.Item3 = new Read10XMerger(r.getReference(), p.Item1);
									}
									p.Item3.add(r, d);
								}
							}
						}
					}
				}
			}
			
//			System.out.println("Assembled reads "+barcodeToPeak.size());
			
			HashMap<GenomicRegion,ArrayList<DefaultAlignedReadsData>> pre = new HashMap<>();
			for (MutablePair<DnaSequence,Integer> pp : regions.keySet()) {
				
				MutableTriple<GenomicRegion, double[], Read10XMerger> p = regions.get(pp);
				
				ImmutableReferenceGenomicRegion<Void> rgr = new ImmutableReferenceGenomicRegion<>(or.getReference(), p.Item1);
				String cellbc = pp.Item1.toString().substring(0,cellBarcodeLength);
				Integer cell = whitelists[pp.Item2].get(cellbc); 
				
				AlignedReadsDataFactory fac = new AlignedReadsDataFactory(nind);
				fac.start();
				fac.newDistinctSequence();
				fac.setMultiplicity(1);
				fac.setCount(cell, 1);
				p.Item3.addVariations(fac);
				
				DefaultAlignedReadsData ard = fac.create();
				pre.computeIfAbsent(rgr.getRegion(), x->new ArrayList<>()).add(ard);
				
			}
//			System.out.println("Return size: "+pre.size());
			return EI.wrap(pre.keySet())
					.map(reg->{
						ArrayList<DefaultAlignedReadsData> l = pre.get(reg);
						if (l.size()==1)
							return new ImmutableReferenceGenomicRegion<>(or.getReference(), reg, l.get(0));
						else {
							AlignedReadsDataFactory fac = new AlignedReadsDataFactory(nind);
							fac.start();
							for (DefaultAlignedReadsData ard : l) 
								fac.add(ard, 0);
							fac.makeDistinct();
							return new ImmutableReferenceGenomicRegion<>(or.getReference(),reg,fac.create());
						}
					});
			
		}




		private static ExtendedIterator<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> demultiplexOld(
				ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> or, 
				ImmutableReferenceGenomicRegion<Void>[] peaks, double mean, double sd, String[] conditions, 
				Genomic g,
				HashMap<String,Integer>[] whitelists) {
			
			Normal dist = new Normal(mean, sd);
			int cellBarcodeLength = whitelists[0].keySet().iterator().next().length();
			int nind = EI.wrap(whitelists).unfold(m->EI.wrap(m.values())).mapToInt(x->x).max()+1;
			
			// collect all reads per cell barcode
			HashMap<Integer,ArrayList<MutablePair<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>,Integer>>> regions = new HashMap<>();
			for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : or.getData()) {
				for (int d=0; d<r.getData().getDistinctSequences(); d++) {
					for (int c=0; c<r.getData().getNumConditions(); c++) {
						for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
							String cell = bc.toString().substring(0,cellBarcodeLength);
							Integer index = whitelists[c].get(cell); 
							if (index!=null) {
								regions.computeIfAbsent(index, x->new ArrayList<>()).add(new MutablePair<>(r, d));
							}
						}
					}
				}
			}
			TrimmedGenomicRegion reuse = new TrimmedGenomicRegion();
			HashMap<GenomicRegion,ArrayList<DefaultAlignedReadsData>> pre = new HashMap<>();
			NumericArray unit = NumericArray.wrap(1);
			for (Integer ind : regions.keySet()) {
				
				double[] majority = new double[peaks.length];
				
				for (MutablePair<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>, Integer> p : regions.get(ind)) {
					for (int i=0; i<majority.length; i++) {
						double w = weight(g,reuse.set(p.Item1.getRegion()),peaks[i],dist);
						majority[i]+=w*p.Item1.getData().getTotalCountForDistinct(p.Item2, ReadCountMode.Weight);
					}
				}
				int m = ArrayUtils.argmax(majority);
				if (majority[m]==0) continue;
				
				GenomicRegion reg = EI.wrap(regions.get(ind)).map(p->extend(g,p.Item1.getRegion(),peaks[m],dist,reuse)).removeNulls().reduce((a,b)->a.union(b));
				ImmutableReferenceGenomicRegion<Void> rgr = new ImmutableReferenceGenomicRegion<>(or.getReference(), reg);
				
				if (rgr.getRegion().getTotalLength()>DefaultAlignedReadsData.MAX_POSITION) {
					reg = rgr.map(new ArrayGenomicRegion(rgr.getRegion().getTotalLength()-DefaultAlignedReadsData.MAX_POSITION,rgr.getRegion().getTotalLength()));
					rgr = new ImmutableReferenceGenomicRegion<>(or.getReference(), reg);
				}
				
				CoverageAlgorithm cov = new CoverageAlgorithm(rgr);
				EI.wrap(regions.get(ind)).forEachRemaining(p->cov.add(p.Item1.getRegion(), unit));
				
				HashMap<Integer,Counter<AlignedReadsVariation>> varsPerPosInRgr = new HashMap<>();
				for (MutablePair<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>,Integer> p : regions.get(ind)) {
					int vc = p.Item1.getData().getVariationCount(p.Item2);
					for (int v=0; v<vc; v++) {
						AlignedReadsVariation vari = p.Item1.getData().getVariation(p.Item2, v);
						if (vari.isDeletion() || vari.isInsertion() || vari.isMismatch()) {
							int pos = p.Item1.map(vari.getPosition());
							if (rgr.getRegion().contains(pos)) { // always except when longer than DefaultAlignedReadsData.MAX_POSITION
								pos = rgr.induce(pos);
								varsPerPosInRgr.computeIfAbsent(pos, x->new Counter<>()).add(vari);
							}
						}
					}
				}
				
				AlignedReadsDataFactory fac = new AlignedReadsDataFactory(nind);
				fac.start();
				fac.newDistinctSequence();
				fac.setMultiplicity(1);
				fac.setCount(ind, 1);
				for (Integer pos : varsPerPosInRgr.keySet()) {
					// majority vote
					Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
					AlignedReadsVariation maxVar = ctr.getMaxElement(0);
					int total = ctr.total()[0];
					if (maxVar!=null && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
						fac.addVariation(maxVar.reposition(pos));
					}
				}
				
				DefaultAlignedReadsData ard = fac.create();
				pre.computeIfAbsent(reg, x->new ArrayList<>()).add(ard);
				
			}
			
			return EI.wrap(pre.keySet())
					.map(reg->{
						ArrayList<DefaultAlignedReadsData> l = pre.get(reg);
						if (l.size()==1)
							return new ImmutableReferenceGenomicRegion<>(or.getReference(), reg, l.get(0));
						else {
							AlignedReadsDataFactory fac = new AlignedReadsDataFactory(nind);
							fac.start();
							for (DefaultAlignedReadsData ard : l) 
								fac.add(ard, 0);
							fac.makeDistinct();
							return new ImmutableReferenceGenomicRegion<>(or.getReference(),reg,fac.create());
						}
					});
			
		}

		private boolean intersects(TrimmedGenomicRegion reg, GenomicRegion read,
				ArrayList<ImmutableReferenceGenomicRegion<Transcript>> ctrans, int maxdist) {
//			return reg.intersects(read.getRegion());
			int tol = maxdist-read.getTotalLength();
			
			if (read.getDistance(reg)<tol) return true;
			
			for (ImmutableReferenceGenomicRegion<Transcript> t : ctrans) {
				if (isCompatible(t.getRegion(), reg)){
					ArrayGenomicRegion iread = t.getRegion().induce(t.getRegion().intersect(read));
					ArrayGenomicRegion ireg = t.getRegion().induce(t.getRegion().intersect(reg));
					if (iread.getDistance(ireg)<tol) return true;
							
				}
			}
			
			return false;
		}

		
		private static GenomicRegion extend(Genomic g, GenomicRegion read,
				ImmutableReferenceGenomicRegion<?> peak, Normal dist, TrimmedGenomicRegion reuse) {
			
			GenomicRegion re = null;
			double w = 0;
			if (peak.isUpstream(read)) {
				w = dist.density(read.getDistance(peak.getRegion()), false);
				if (peak.getReference().isMinus())
					re = new ArrayGenomicRegion(peak.getRegion().getStart(),read.getStart()).union(read);
				else
					re = new ArrayGenomicRegion(read.getEnd(),peak.getRegion().getStart()).union(read);
			}
			if (read.contains(peak.getRegion()))
				throw new RuntimeException(read+" "+peak);

			for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(peak).chain(g.getTranscripts().ei(peak.toMutable().toOppositeStrand()))
					.filter(t->isCompatible(t.getRegion(),reuse.set(read))).loop()) {
				GenomicRegion mr = induce(t, reuse.set(read));
				int mp = t.induce(peak.getRegion().getStart());
				if (mp>mr.getStop()) {
					double nw = dist.density(mp-mr.getStop(), false);
					if (nw>w) {
						w = nw;
						re = t.map(new ArrayGenomicRegion(mr.getEnd(),mp).union(mr));
					}
				}
				
			}

			return re;
			
		}

		private static double weight(Genomic g, TrimmedGenomicRegion read,
				ImmutableReferenceGenomicRegion<?> peak, Normal dist) {
			double w = 0;
			if (peak.isUpstream(read.parent))
				w = dist.density(read.parent.getDistance(peak.getRegion()), false);
			if (read.parent.contains(peak.getRegion()))
				return 0;
			
			double max = g.getTranscripts().ei(peak).chain(g.getTranscripts().ei(peak.toMutable().toOppositeStrand()))
				.filter(t->isCompatible(t.getRegion(), read))
				.mapToDouble(t->{
					GenomicRegion mr = induce(t, read);
					int mp = t.induce(peak.getRegion().getStart());
					if (mp>mr.getStop()) return dist.density(mp-mr.getStop(), false);
					return 0;
				}).max();
			
			return Math.max(max, w);
		}
		
	}
	
	public static class ProfileWalker {
		private int dist;
		private int pos;
		private double value = 0;
		private int lastpos = -1;
		private double lastvalue = 0;
		private int direction;
		private boolean withread;
		private double[] sm;
		private double[] sm2;
		private boolean pathHasRead;
		
		private double bestDiffVal = Double.NEGATIVE_INFINITY;
		private int bestDiffPos = -1;
		
		@Override
		public String toString() {
			return pos+" "+direction+" "+withread;
		}
		
		public ProfileWalker(int pos, int direction, double[] sm, double[] sm2) {
			this.pos = pos;
			this.direction = direction;
			value = sm[pos];
			this.sm = sm;
			this.sm2 = sm2;
			if (sm2!=null && sm2[pos]>0) withread=true;
			pathHasRead = sm2!=null && sm2[pos]>0;
		}

		public ProfileWalker(int dist, int pos, double value, int lastpos, double lastvalue, int direction, boolean withread, double[] sm,
						double[] sm2, double bestDiffVal, int bestDiffPos) {
			this.dist = dist;
			this.pos = pos;
			this.value = value;
			this.lastpos = lastpos;
			this.lastvalue = lastvalue;
			this.direction = direction;
			this.withread = withread;
			this.sm = sm;
			this.sm2 = sm2;
			this.bestDiffVal = bestDiffVal;
			this.bestDiffPos = bestDiffPos;
		}
		public boolean moveAndCheckIntrons(SpliceGraph introns, Collection<ProfileWalker> add, int maxo) {
			
			lastpos = pos;
			lastvalue = value;
			pos+=direction;
			if (pos<0 || pos>=sm.length) return false;
			
			value = sm[pos];
			sm[pos] = 0;
			
			pathHasRead |= sm2!=null && sm2[pos]>0;
			int before = add.size();
			
			introns.getIntrons().forEachIntervalIntersecting(pos,pos,intron->{
				if (lastpos<intron.getStart()) {
					ProfileWalker w = new ProfileWalker(dist+1,intron.getEnd(),sm[intron.getEnd()],lastpos,lastvalue,1,withread,sm,sm2,bestDiffVal,bestDiffPos);
					if (w.isValid(maxo)) {
						add.add(w);
						sm[w.pos]=0;
					}
					w = new ProfileWalker(dist+1,intron.getStop(),sm[intron.getStop()],lastpos,lastvalue,-1,false,sm,null,Double.NEGATIVE_INFINITY,-1);
					if (w.isValid(0)) {
						add.add(w);
						sm[w.pos]=0;
					}
				}
				else if (lastpos>intron.getStop()) {
					ProfileWalker w = new ProfileWalker(dist+1,intron.getStart(),sm[intron.getStart()],lastpos,lastvalue,1,false,sm,null,Double.NEGATIVE_INFINITY,-1);
					if (w.isValid(0)) {
						add.add(w);
						sm[w.pos]=0;
					}
					w = new ProfileWalker(dist+1,intron.getStart()-1,sm[intron.getStart()-1],lastpos,lastvalue,-1,false,sm,null,Double.NEGATIVE_INFINITY,-1);
					if (w.isValid(maxo)) {
						add.add(w);
						sm[w.pos]=0;
					}
				}
			});
			
			if (add.size()>before)
				pathHasRead = false;

			dist++;
			withread &= sm2!=null && sm2[pos]>0; 
			
			if (dist==maxo || (withread && !inRange(maxo))) {
				double curdiff = lastvalue-value;
				if (curdiff>0.5*bestDiffVal) {
					bestDiffPos = lastpos;
					bestDiffVal = curdiff;
				}
			}
			
			return isValid(maxo);
			
		}
		
		private boolean hasPathRead() {
			return pathHasRead;
		}
		
		private boolean isValid(int maxo) {
			return inRange(maxo) || (!increased() && !isNull());
		}

		private boolean inRange(int maxo) {
			return dist<=maxo;
		}
		
		private boolean isNull() {
			return value==0;
		}

		private boolean increased() {
			return value>lastvalue;
		}
		
	}
	
	private static NumericArray unit = NumericArray.wrap(1);
	public static class Read10XMerger {

		
		HashMap<Integer,Counter<AlignedReadsVariation>> varsPerPosInRgr;
		CoverageAlgorithm cov;
		
		public Read10XMerger(ReferenceSequence ref, GenomicRegion region) {
			cov = new CoverageAlgorithm(ref,region);
			varsPerPosInRgr = new HashMap<>();
		}
		
		public void addVariations(AlignedReadsDataFactory fac) {
			for (Integer pos : varsPerPosInRgr.keySet()) {
				// majority vote
				Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
				AlignedReadsVariation maxVar = ctr.getMaxElement(0);
				int total = ctr.total()[0];
				if (maxVar!=null && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
					fac.addVariation(maxVar.reposition(pos));
				}
			}
		}

		public void add(ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r, int d) {
			cov.add(r.getRegion(), unit);
			int vc = r.getData().getVariationCount(d);
			for (int v=0; v<vc; v++) {
				AlignedReadsVariation vari = r.getData().getVariation(d, v);
				if (vari.isDeletion() || vari.isInsertion() || vari.isMismatch()) {
					int pos = r.map(vari.getPosition());
					if (cov.getParentRegion().getRegion().contains(pos)) { // always except when longer than DefaultAlignedReadsData.MAX_POSITION
						pos = cov.getParentRegion().induce(pos);
						varsPerPosInRgr.computeIfAbsent(pos, x->new Counter<>()).add(vari);
					}
				}
			}
		}
		
	}

	public static class Read10XEmObject {
		ImmutableReferenceGenomicRegion<Void>[] peaks;
		BitVector necessaryPeak;
		int[] usePeaks;
		double[] peaksum;
		
		public Read10XEmObject(ImmutableReferenceGenomicRegion<Void>[] peaks) {
			this.peaks = peaks;
			necessaryPeak = new BitVector(peaks.length);
			peaksum = new double[peaks.length];
			usePeaks = EI.seq(0, peaks.length).toIntArray();
		}

		public <T> HashMap<T, ImmutableReferenceGenomicRegion<Void>> getMap(
				HashMap<T, MutableTriple<GenomicRegion, double[],Read10XMerger>> regions) {
			HashMap<T, ImmutableReferenceGenomicRegion<Void>> re = new HashMap<T, ImmutableReferenceGenomicRegion<Void>>();
			for (T bc : regions.keySet()) {
				
				MutableTriple<GenomicRegion, double[],Read10XMerger> rdp = regions.get(bc);
				double max = 0;
				int maxi = -1;
				for (int i : usePeaks) {
					double h = rdp.Item2[i]*peaksum[i];
					if (h>max) {
						max = h;
						maxi = i;
					}
				}
				if (maxi==-1) 
					throw new RuntimeException();
				re.put(bc, peaks[maxi]);
			}
			return re;
		}

		public boolean add(Genomic g, Normal dist,MutableTriple<GenomicRegion,double[],Read10XMerger> rdp, TrimmedGenomicRegion reuse) {
			double s = 0;
			int necIndex = -1;
			for (int i : usePeaks) {
				double w = Read10XClusterProgram.weight(g,reuse.set(rdp.Item1),peaks[i],dist);
				if (w>0) {
					rdp.Item2[i] = w;
					s+=rdp.Item2[i];
					if (necIndex==-1) necIndex=i;
					else necIndex=-2;
				}
			}
			if (necIndex==-1) 
				return false;
			
			for (int i : usePeaks) { 
				rdp.Item2[i]=rdp.Item2[i]/s;
				peaksum[i] += rdp.Item2[i];
			}
			if (necIndex>=0) necessaryPeak.putQuick(necIndex,true);
			return true;
		}

		public boolean checkUnnecessary(Collection<MutableTriple<GenomicRegion,double[],Read10XMerger>> regions, double thres) {
			boolean changed = true;
			boolean anychanged = false;
			while (changed) {
				BitVector hereNecessaryPeak = necessaryPeak.copy();
				for (MutableTriple<GenomicRegion,double[],Read10XMerger> rdp : regions) {
					double max = 0;
					for (int i : usePeaks) 
						max = Math.max(max,rdp.Item2[i]*peaksum[i]);
					for (int i : usePeaks)
						if (rdp.Item2[i]*peaksum[i]>=max*thres)
							hereNecessaryPeak.putQuick(i,true);
				}
				
				changed = hereNecessaryPeak.cardinality()<usePeaks.length;
				if (changed) {
					// only remove one at a time! select the weakest one
					int weakest = -1;
					for (int i : usePeaks) 
						if (!hereNecessaryPeak.getQuick(i) && (weakest==-1 || peaksum[i]<peaksum[weakest]))
							weakest=i;
					
					usePeaks = ArrayUtils.removeIndexFromArray(usePeaks, ArrayUtils.find(usePeaks,weakest));
					peaksum[weakest] = 0;
					
					anychanged = true;
				}
			}
			
			return anychanged;
		}

		public void iterate(Collection<MutableTriple<GenomicRegion,double[],Read10XMerger>> regions, int rounds) {
			ArrayUtils.normalize(peaksum);
			
			double[] npeaksum = new double[peaksum.length];
			for (int iter=0; iter<rounds; iter++) {
				for (MutableTriple<GenomicRegion,double[],Read10XMerger> rdp : regions) {
					double s = 0;
					for (int i : usePeaks) 
						s+=rdp.Item2[i]*peaksum[i];
					for (int i : usePeaks)  
						npeaksum[i]+=rdp.Item2[i]*peaksum[i]/s;
				}
				System.arraycopy(npeaksum, 0, peaksum, 0, peaksum.length);
				ArrayUtils.normalize(peaksum);
				Arrays.fill(npeaksum, 0);
			}			
		}
	}
	
	public static class Read10XInfer3pProgram extends GediProgram {

		public Read10XInfer3pProgram(Read10XParameterSet params) {
			addInput(params.prefix);
			addInput(params.reads);
			addInput(params.genomic);
			addInput(params.nthreads);
			
			
			addOutput(params.infer3pparam);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			GenomicRegionStorage<BarcodedAlignedReadsData> reads = getParameter(1);
			Genomic genomic = getParameter(2);
			int nthreads = getIntParameter(3);

			
			String[] signals = {"AATAAA","ATTAAA","AGTAAA","TATAAA"};
			Trie<Integer> strie = new Trie<>();
			for (int i=0; i<signals.length-1; i++)
				strie.put(signals[i], i);
			
			int[] signalrange = {12,30}; // most of the ccds transcripts have the start of the aataaa in there! 
			int minDistOther = 1000;
			
			int[] count = new int[6];
			
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> niceTranscripts = genomic.getTranscripts().ei()
				.progress(context.getProgress(), (int)genomic.getTranscripts().size(), t->t.getData().getTranscriptId())
				.sideEffect(x->count[0]++)
				.filter(t->t.getData().isCoding())
				.sideEffect(x->count[1]++)
				.filter(t->{
					int t3p = GenomicRegionPosition.ThreePrime.position(t);
					int t3p1 = GenomicRegionPosition.ThreePrime.position(t,1);
					return genomic.getTranscripts()
						.ei(new ImmutableReferenceGenomicRegion<>(t.getReference(), new ArrayGenomicRegion(Math.min(t3p1,t3p),Math.max(t3p1,t3p))))
						.filter(o->!o.equals(t))
						.filter(o->Math.abs(GenomicRegionPosition.ThreePrime.position(o)-t3p)<minDistOther)
						.count()==0;
				})
				.sideEffect(x->count[2]++)
				.map(t->new MutableTriple<>(t,genomic.getSequence(t.getData().get3Utr(t)).toString(),genomic.getSequence(t.getDownstream(50)).toString()))
				.filter(tr->{
					String expectPas = tr.Item2.substring(Math.max(0, tr.Item2.length()-signalrange[1]));
					int paspos = expectPas.indexOf(signals[0]);
					return paspos>=0 && paspos<signalrange[1]-signalrange[0];
				})
				.sideEffect(x->count[3]++)
				.filter(tr->tr.Item2.length()>signalrange[1]?
					strie.iterateAhoCorasick(tr.Item2.substring(0,tr.Item2.length()-signalrange[1])).count()==0
					:true
				)
				.sideEffect(x->count[4]++)
				.filter(tr->
					getFirstPolyAStretch(tr.Item2+tr.Item3).size()==0
				)
				.sideEffect(x->count[5]++)
				.map(tr->tr.Item1)
				.list();
			
			context.logf("Initial transcripts:  %d",count[0]);
			context.logf("Removed non-coding:   %d",count[1]);
			context.logf("Removed unclear 3':   %d",count[2]);
			context.logf("Removed no PAS:       %d",count[3]);
			context.logf("Removed another PAS:  %d",count[4]);
			context.logf("Removed pA stretches: %d",count[5]);
			
			Collections.shuffle(niceTranscripts);
			
			WeightedMeanVarianceOnline mvo = EI.wrap(niceTranscripts)
				.parallelized(nthreads, 8, ei->ei.map(t->{
					TrimmedGenomicRegion tgr = new TrimmedGenomicRegion();
					BarcodeCounter bco = new BarcodeCounter();
					for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : reads.ei(t).filter(r->isCompatible(t.getRegion(), tgr.set(r.getRegion()))).loop()) {
						for (int c=0; c<read.getData().getNumConditions(); c++)
							for (int d=0; d<read.getData().getDistinctSequences(); d++)
								bco.addBarcodes(read.getData().getBarcodes(d, c));
					}
					HashMap<DnaSequence, MutableInteger> bcoo = bco.get();
					WeightedMeanVarianceOnline tmvo = new WeightedMeanVarianceOnline();
					for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : reads.ei(t).filter(r->isCompatible(t.getRegion(), tgr.set(r.getRegion()))).loop()) {
						int dist = t.getRegion().getTotalLength()-induce(t,tgr.set(read.getRegion())).getStart();
						if (dist<=450) {
							double w = 0;
							for (int c=0; c<read.getData().getNumConditions(); c++)
								for (int d=0; d<read.getData().getDistinctSequences(); d++)
									for (DnaSequence bc : read.getData().getBarcodes(d, c)) 
										w+=1.0/bcoo.get(bc).N;
							tmvo.add(dist, w);
						}
					}
					return new MutablePair<>(tmvo,t);
				}))
				.progress(context.getProgress(),(int)niceTranscripts.size(), t->t.Item2.toLocationString())
				.reduce(new WeightedMeanVarianceOnline(),(a,b)->b.add(a.Item1));

			FileUtils.writeAllLines(new String[] {"Name\tValue","Mean\t"+mvo.getMean(),"Sd\t"+mvo.getStandardDeviation(),""}, getOutputFile(0));

			return null;
		}
	}

	
	private static boolean isCompatible(GenomicRegion t,TrimmedGenomicRegion read) {
//		GenomicRegion rs = read.extendBack(-3).extendFront(-3);
		return //t.containsUnspliced(rs);
				t.intersect(read).getNumParts()==read.getNumParts() 
				&& t.isIntronConsistent(read);
	}


//	private static GenomicRegion induce(ReferenceGenomicRegion<?> t,ReferenceGenomicRegion<?> r) {
//		return induce(t,r.getRegion());
//	}
	private static GenomicRegion induce(ReferenceGenomicRegion<?> t,TrimmedGenomicRegion r) {
		return t.induce(r).extendFront(3).extendBack(3).intersect(new ArrayGenomicRegion(0,t.getRegion().getTotalLength()));
	}

	private static class TrimmedGenomicRegion implements GenomicRegion {
		private GenomicRegion parent;
		
		public TrimmedGenomicRegion set(GenomicRegion reg) {
			this.parent = reg;
			return this;
		}
		
		@Override
		public String toString() {
			return toString2();
		}
		
		@Override
		public int getStart(int part) {
			if (!trimStart()) 
				return part==0?parent.getStart()+3:parent.getStart(part);
			return parent.getStart(part+1);
		}
		
		@Override
		public int getEnd(int part) {
			if (part==getNumParts()-1)
				return trimEnd()?parent.getEnd(parent.getNumParts()-2):parent.getEnd()-3;
				
			return parent.getEnd(trimStart()?part+1:part);
		}
		
		private final boolean trimStart() {
			return parent.getLength(0)<=3;
		}
		
		private final boolean trimEnd() {
			return parent.getLength(parent.getNumParts()-1)<=3;
		}
		
		@Override
		public int getNumParts() {
			return parent.getNumParts()-(trimStart()?1:0)-(trimEnd()?1:0);
		}
		
	}
	
	private static class BarcodeCounter {
		private HashMap<DnaSequence,MutableInteger> counter = new HashMap<DnaSequence,MutableInteger>();
		
		public void addBarcodes(DnaSequence[] dna) {
			for (DnaSequence d : dna) {
				counter.computeIfAbsent(d, x->new MutableInteger()).N++;
			}
		}
		
		public HashMap<DnaSequence, MutableInteger> get() {
			return counter;
		}
		
//		public HashMap<DnaSequence,Double> estimate() {
//			if (counter.elements().size()==0) return new HashMap<>();
//			
//			HashMap<DnaSequence, Double> re = new HashMap<>();
//			// deal with 1-hamming distance
//			SimpleDirectedGraph<DnaSequence> mg = new MismatchGraphBuilder<>(counter.elements().toArray(new DnaSequence[0])).build();
//			for (SimpleDirectedGraph<DnaSequence> c : mg.getWeaklyConnectedComponents()) {
//				double sum = EI.wrap(c.getSources()).mapToDouble(d->counter.get(d)[0]).sum();
//				for (DnaSequence d : c.getSources())
//					re.put(d, sum);
//			}
//			
//			for (DnaSequence d : counter.elements())
//				if (!re.containsKey(d))
//					re.put(d, counter.get(d)[0]+0.0);
//			
//			return re;
//		}
		
	}

	
	
	public static class Read10XParameterSet extends GediParameterSet {

		public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
		public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "Output prefix", true, new StringParameterType());
		public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"g", "Genomic name", true, new GenomicParameterType());
		public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "The mapped reads from the ribo-seq experiment.", false, new StorageParameterType<AlignedReadsData>());
		public GediParameter<Integer> maxdist = new GediParameter<Integer>(this,"maxdist", "The maximal distance between two read 5' ends", false, new IntParameterType(), 200);
		public GediParameter<String> whitelist = new GediParameter<String>(this,"cells", "File containing the cell barcodes to use", false, new StringParameterType());
		public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Produce plots", false, new BooleanParameterType());

		public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "Only search against chromosome 22", false, new BooleanParameterType());
		
		
		public GediParameter<File> clustercit = new GediParameter<File>(this,"${prefix}.cluster.cit", "Clusters for the viewer", false, new FileParameterType());
		public GediParameter<File> infer3ptsv = new GediParameter<File>(this,"${prefix}.3p.tsv", "Table of 3' end information.", false, new FileParameterType());
		public GediParameter<File> infer3pparam = new GediParameter<File>(this,"${prefix}.3p.parameter", "3' parameter file.", false, new FileParameterType());
		
	}

	
}
