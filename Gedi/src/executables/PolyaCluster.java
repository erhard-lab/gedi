package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.NavigableMap;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsSoftclip;
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
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.MemoryIntegerArray;
import gedi.util.datastructure.array.SparseMemoryCountArray;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.datastructure.tree.Trie.AhoCorasickResult;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.datastructure.tree.redblacktree.UnitInterval;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;


public class PolyaCluster {
	
// for Lexogen 3' (FWD) data!




	public static void main(String[] args) throws IOException {
		
		PolyaParameterSet params = new PolyaParameterSet();
		GediProgram pipeline = GediProgram.create("Polya",
//				new PolyaClusterPositionInferenceProgram(params),
				new PolyaClusterProgram(params),
				new PolyaAnnotateClusterProgram(params),
				new PolyaAnalyzeClusterProgram(params)
				);
		GediProgram.run(pipeline, null, new CommandLineHandler("Polya","Polya analyzes 3' seq data.",args));

	}
	
	public enum PeakCategory {
		Exonic,Spliced,ExonIntronic,Intronic,Flank,Intergenic
	}

	public static class PolyaAnalyzeClusterProgram extends GediProgram {

		public PolyaAnalyzeClusterProgram(PolyaParameterSet params) {
			addInput(params.prefix);
			addInput(params.genomic);
			addInput(params.annotation);
			
			addOutput(params.pairwise);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			Genomic g = getParameter(1);
			File tab = getParameter(2);
			
			
			HeaderLine h = new HeaderLine();
			ToIntFunction<String[]> toTotal = a->Integer.parseInt(a[h.get("Total")]);
			Predicate<String[]> isNormal = a->a[h.get("ThreePrime")].equals("1");
			Predicate<String[]> isInternal = a->!a[h.get("ThreePrime")].equals("1") && !a[h.get("Signal")].equals("-") && a[h.get("Priming")].equals("0") && a[h.get("Antisense")].equals("0");
			Comparator<String[]> gcomp = (a,b)->a[h.get("Gene")].compareTo(b[h.get("Gene")]);
			
			LineWriter out = null;
			for (String[][] gene : EI.lines(tab.getPath()).header(h).split("\t")
										.sort(gcomp)
										.progress(context.getProgress(), -1, d->d[0])
										.multiplex(gcomp, String[].class).loop()) {
				
				if (out==null) {
					out = getOutputWriter(0);
					out.writeLine("Symbol\tLocation\t"+EI.wrap(h.getFields()).map(f->f+".A").concat("\t")+"\t"+EI.wrap(h.getFields()).map(f->f+".B").concat("\t"));
				}
				String symbol = g.getGeneTable("symbol").apply(gene[0][h.get("Gene")]);
				
				if (!"null".equals(gene[0][h.get("Gene")])) {
					Arrays.sort(gene, (a,b)->{
						int re=(isNormal.test(a)?0:1)-(isNormal.test(b)?0:1);
						if (re!=0) return re;
						return Integer.compare(toTotal.applyAsInt(b),toTotal.applyAsInt(a));
					});
					
					if (isNormal.test(gene[0]))
						for (int i=1; i<gene.length; i++) {
							MutableReferenceGenomicRegion<Object> loc = ImmutableReferenceGenomicRegion.parse(gene[0][h.get("Location")]).toMutable();
							loc.setRegion(loc.getRegion().extendBack(150).extendFront(150).union(ImmutableReferenceGenomicRegion.parse(gene[i][h.get("Location")]).getRegion().extendBack(150).extendFront(150)));
							out.writef("%s\t%s\t%s\t%s\n", symbol,loc.toLocationString(),StringUtils.concat("\t", gene[0]),StringUtils.concat("\t", gene[i]));
						}

				}
				
			}
			if (out!=null) 
				out.close();
				
			return null;
		}
	}
	
	public static class PolyaAnnotateClusterProgram extends GediProgram {

		public PolyaAnnotateClusterProgram(PolyaParameterSet params) {
			addInput(params.prefix);
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.clustercit);
			addInput(params.reads);
			
			addOutput(params.annotation);
			
			
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			Genomic g = getParameter(1);
			int nthreads = getIntParameter(2);
			File cl = getParameter(3);
			GenomicRegionStorage<DefaultAlignedReadsData> reads = getParameter(4);
			
			int downstreamAnnotated = 120;
			int tol = 50;
			
			String[] signals = {"AATAAA","ATTAAA","AGTAAA","TATAAA","-"};
			Trie<Integer> strie = new Trie<>();
			for (int i=0; i<signals.length-1; i++)
				strie.put(signals[i], i);
		
			String[] conditions = reads.getMetaDataConditions();
			int gene_tolerance = 10_000;
			
			CenteredDiskIntervalTreeStorage<SparseMemoryCountArray> clusters = new CenteredDiskIntervalTreeStorage<>(cl.getPath());
			
			MemoryIntervalTreeStorage<String> extendedGenes = new MemoryIntervalTreeStorage<>(String.class);
			extendedGenes.fill(g.getGenes().ei().map(ge->ge.toMutable().alterRegion(reg->reg.extendBack(gene_tolerance).extendFront(gene_tolerance)).toImmutable()));
			
			
			LineWriter out = getOutputWriter(0);
			out.writef("Location\tGene\tCategory\tThreePrime\tAntisense\tSignal\tSignal.dist\tPriming\tTotal\n");
			
			for (ImmutableReferenceGenomicRegion<SparseMemoryCountArray> p : clusters.ei().progress(context.getProgress(), (int)clusters.size(), d->d.toLocationString()).loop()) {
				
				PeakCategory cat = PeakCategory.Intergenic;
				String gene = null;
				boolean threeprime = false;
				boolean antisense = false;

				for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(p.getReference(),p.getRegion().extendBack(tol).extendFront(tol)).loop()) {
					if (t.getRegion().extendBack(tol).extendFront(tol).containsUnspliced(p.getRegion())) {
						cat = p.getRegion().getNumParts()==1?PeakCategory.Exonic:PeakCategory.Spliced;
						gene = t.getData().getGeneId();
						if (t.induceMaybeOutside(p.getReference().isMinus()?p.getRegion().getStart():p.getRegion().getStop())+downstreamAnnotated*2>t.getRegion().getTotalLength())
							threeprime = true;
					}
				}
				if (cat==PeakCategory.Intergenic) {
					for (ImmutableReferenceGenomicRegion<String> gg : g.getGenes().ei(p).loop()) {
						if (gg.getRegion().contains(p.getRegion())) {
							boolean exin = g.getTranscripts().ei(gg).filter(t->t.getRegion().intersects(p.getRegion())).count()>0;
							PeakCategory ccat = exin?PeakCategory.ExonIntronic:PeakCategory.Intronic;
							if (cat.ordinal()>ccat.ordinal()) cat = ccat;
							gene = gg.getData();
						}
					}
				}
				if (cat==PeakCategory.Intergenic) {
					if (extendedGenes.ei(p).count()>0) {
						cat = PeakCategory.Flank;
						gene = extendedGenes.ei(p).first().getData();
					}
				}
				
				if (cat==PeakCategory.Intergenic) {
					for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(p.toMutable().toOppositeStrand()).loop()) {
						if (t.getRegion().containsUnspliced(p.getRegion())) {
							cat = t.getRegion().getNumParts()==1?PeakCategory.Exonic:PeakCategory.Spliced;
							gene = t.getData().getGeneId();
							antisense = true;
						}
					}
				}
				if (cat==PeakCategory.Intergenic) {
					for (ImmutableReferenceGenomicRegion<String> gg : g.getGenes().ei(p.toMutable().toOppositeStrand()).loop()) {
						if (gg.getRegion().contains(p.getRegion())) {
							boolean exin = g.getTranscripts().ei(gg).filter(t->t.getRegion().intersects(p.getRegion())).count()>0;
							PeakCategory ccat = exin?PeakCategory.ExonIntronic:PeakCategory.Intronic;
							if (cat.ordinal()>ccat.ordinal()) cat = ccat;
							gene = gg.getData();
							antisense = true;
						}
					}
				}
				if (cat==PeakCategory.Intergenic) {
					if (extendedGenes.ei(p.toMutable().toOppositeStrand()).count()>0) {
						cat = PeakCategory.Flank;
						gene = extendedGenes.ei(p.toMutable().toOppositeStrand()).first().getData();
						antisense = true;
					}
				}
				
				
				int fs = signals.length-1;
				int dist = 0;
				
				CharSequence upstream = g.getSequenceSave(p.getReference(),p.mapMaybeOutSide(new ArrayGenomicRegion(-30,-5)));
				CharSequence sseq = StringUtils.reverse(g.getSequenceSave(p.getReference(),p.mapMaybeOutSide(new ArrayGenomicRegion(-15,1)))).toString();
				
				int as = 0; int gs = 0; int oths = 0;
				int max = 0;
				for (int i=0; i<sseq.length(); i++) {
					if (sseq.charAt(i)=='A') as++;
					else if (sseq.charAt(i)=='G') gs++;
					else oths++;
					if (oths>1) {
						break;
					}
					if (as>gs) max=as+gs;
				}
				
				boolean stretch = max>=6;
				
				
				for (AhoCorasickResult<Integer> hit : strie.iterateAhoCorasick(upstream).loop()) {
					if (hit.getValue()<fs) {
						fs = hit.getValue();
						dist = 30-hit.getStart();
					}
				}
				
				out.writef("%s\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t%.0f\n", 
						p.toLocationString(),
						gene,
						cat,
						threeprime?1:0,
						antisense?1:0,
						signals[fs],
						dist,
						stretch?1:0,
						p.getData().sum());
				
			}
			
			out.close();
			
			if (conditions.length>8) {
				File dir = new File(prefix+".quant");
				dir.mkdirs();
				if (!dir.isDirectory()) throw new IOException("Could not create directory "+dir);
				
				EI.wrap(conditions).print(new File(dir,"barcodes.tsv.gz").getPath());
				out = new LineOrientedFile(dir, "features.tsv.gz").write();
				int nonzeros = 0;
				for (ImmutableReferenceGenomicRegion<SparseMemoryCountArray> p : clusters.ei().progress(context.getProgress(), (int)clusters.size(), d->d.toLocationString()).loop()) {
					out.writef("%s\t%s\tGene Expression\n", g,p.toLocationString());
					nonzeros+=p.getData().getNonZeroCount();
				}
					
				out.close();
				
				LineWriter writer = new LineOrientedFile(dir,"matrix.mtx.gz").write();
				writer.writeLine("%%MatrixMarket matrix coordinate real general");
				writer.writeLine("%metadata_json: {\"format_version\": 2, \"software_version\": \"3.0.2\"}");
				writer.writef("%d %d %d\n", clusters.size(),conditions.length,nonzeros);
				
				MutableInteger ge = new MutableInteger(1);
				for (ImmutableReferenceGenomicRegion<SparseMemoryCountArray> p : clusters.ei().progress(context.getProgress(), (int)clusters.size(), d->d.toLocationString()).loop()) {
					p.getData().forEachNonZero((idx,val)->writer.writef2("%d %d %d\n", ge.N, idx+1,p.getData().getInt(idx)));
					ge.N++;
				}
				writer.close();
				
			}
			else {
				LineWriter writer = new LineOrientedFile(prefix+".quant.tsv.gz").write();
				writer.writeLine("Location+\t"+EI.wrap(conditions).concat("\t"));
				for (ImmutableReferenceGenomicRegion<SparseMemoryCountArray> p : clusters.ei().progress(context.getProgress(), (int)clusters.size(), d->d.toLocationString()).loop()) {
					writer.writef("%s\t%s\n",p.toLocationString(),p.getData().formatArray("\t"));
				}
				writer.close();
			}
			
			return null;
		}
	}
	
	private static class PolyaStat {
		MemoryIntegerArray inpolya;
		MemoryIntegerArray nonpolya;
		MemoryIntervalTreeStorage<MutableInteger> nonRegions;
		
		public PolyaStat(int cond) {
			inpolya = new MemoryIntegerArray(cond);
			nonpolya = new MemoryIntegerArray(cond);
			nonRegions = new MemoryIntervalTreeStorage<>(MutableInteger.class);
		}

		public PolyaStat add(PolyaStat ostat) {
			inpolya.add(ostat.inpolya);
			nonpolya.add(ostat.nonpolya);
			nonRegions.fill(ostat.nonRegions.ei());
			
			return this;
		}
	}
	
	public static class PolyaClusterProgram extends GediProgram {



		public PolyaClusterProgram(PolyaParameterSet params) {
			addInput(params.prefix);
			addInput(params.reads);
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.maxFrag);
			addInput(params.minPrimerLen);
			addInput(params.minPrimerFrac);
			addInput(params.mergeDist);
			addInput(params.test);
			addInput(params.minReads);
			addInput(params.minConditions);
			
			addOutput(params.clustercit);
			addOutput(params.readStat);
//			addOutput(params.nonCit);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			GenomicRegionStorage<DefaultAlignedReadsData> reads = getParameter(1);
			Genomic genomic = getParameter(2);
			int nthreads = getIntParameter(3);
			int maxFrag = getIntParameter(4);
			int minPrimerLen  = getIntParameter(5);
			double minPrimerFrac = getDoubleParameter(6);
			int mergeDist = getIntParameter(7);
			boolean test = getBooleanParameter(8);
			int minReads = getIntParameter(9);
			int minCondp = getIntParameter(10);
			
			genomic.getTranscriptTable("source");
			int numCond = reads.getRandomRecord().getNumConditions();
			
			int minCond = minCondp<0?(int)Math.ceil(numCond*(-minCondp)/100):minCondp;
			
			context.getLog().info("Identifying clusters...");
			
			
			TreeSet<ReferenceSequence> refs = new TreeSet<ReferenceSequence>(reads.getReferenceSequences());
			refs.removeIf(r->!genomic.getSequenceNames().contains(r.getName()));
			
			CenteredDiskIntervalTreeStorage<SparseMemoryCountArray> out = new CenteredDiskIntervalTreeStorage<>(getOutputFile(0).getPath(),SparseMemoryCountArray.class);
			
			Predicate<ImmutableReferenceGenomicRegion<SparseMemoryCountArray>> filter = rgr->{
				MutableInteger c = new MutableInteger();
				rgr.getData().forEachNonZero((idx,val)->{
					if (val>=minReads)
						c.N++;
				});
				return c.N>=minCond;
			};
			
			PolyaStat allstat = new PolyaStat(numCond);
			int block = 16;
			IterateIntoSink<ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>>> sink = new IterateIntoSink<>(
					ei->out.fill(ei.parallelized(nthreads, block, 
										()->new PolyaStat(numCond),
										(xei,stat)->xei.unfold(x->findPolya(genomic, x, maxFrag, minPrimerLen, minPrimerFrac, mergeDist, stat)	)
										).executeStatesWhenFinished(statei->statei.forEachRemaining(ostat->allstat.add(ostat)))
										.removeNulls()
										.filter(filter)
										.progress(context.getProgress(), -1, x->"Detecting Poly-A sites "+x.toLocationString())
										,context.getProgress()),
					nthreads*block
					);
			
			
			IntervalTree<GenomicRegion,ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>> openClusters = null;
			TreeMap<Integer,GenomicRegion> openClustersByStop = null;
			
			MemoryIntervalTreeStorage<Transcript> transcripts = genomic.getTranscripts();
			
			for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read : (test?reads.ei("22+:20586947-20588174").sort():reads.ei()) //reads.ei("10+")   
										.progress(context.getProgress(), (int)reads.size(), r->"Collecting reads "+r.toLocationString())
										.filter(r->r.getData().getMultiplicity(0)<=1)
										.checkOrder((a,b)->a.compareTo(b))
										.loop()) {
				
				if (read.getData().getNumParts(read,0)>1 && transcripts.ei(read).chain(transcripts.ei(read.toMutable().toOppositeStrand()))
							.filter(t->read.getData().isConsistentlyContained(read, t, 0))
							.count()==0)
					continue;
				
				ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> re;
				if (openClusters==null) {
					// first cluster
					openClusters = new IntervalTree<GenomicRegion, ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>>(read.getReference());
					openClusters.put(read.getRegion(),re = new ArrayList<>());
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				}
				else if (!openClusters.getReference().equals(read.getReference())) {
					// first cluster on this chromosome
					
					for (GenomicRegion clreg : openClusters.keySet()) 
						sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, openClusters.get(clreg)));

					openClusters = new IntervalTree<GenomicRegion, ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>>(read.getReference());
					openClusters.put(read.getRegion(),re = new ArrayList<>());
					openClustersByStop = new TreeMap<>();
					openClustersByStop.put(read.getRegion().getStop(), read.getRegion());
				} else {
					ArrayList<GenomicRegion> overlaps = openClusters
						.keys(read.getRegion().getStart(), read.getRegion().getStop())
						.filter(reg->reg.intersects(read.getRegion())).list();
					
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
			
			if (openClusters!=null)
				for (GenomicRegion clreg : openClusters.keySet()) 
					sink.accept(new ImmutableReferenceGenomicRegion<>(openClusters.getReference(), clreg, openClusters.get(clreg)));

			
			sink.finish();
			
			LineWriter rstat = getOutputWriter(1);
			rstat.writeLine("Condition\tInside\tOutside");
			for (int i=0; i<numCond; i++)
				rstat.writef("%s\t%d\t%d\n",reads.getMetaDataConditions()[i], allstat.inpolya.getInt(i), allstat.nonpolya.getInt(i));
			rstat.close();

			//context.logf("Writing read clusters non associated with polyA.");
			//new CenteredDiskIntervalTreeStorage<>(getOutputFile(2).getPath(),MutableInteger.class).fill(allstat.nonRegions,context.getProgress());
			
			
			return null;
		}
		
		private static AlignedReadsSoftclip getSoftclip(DefaultAlignedReadsData ard, int distinct) {
			for (int v=0; v<ard.getVariationCount(distinct); v++)
				if (ard.isSoftclip(distinct, v) && !ard.isSoftclip5p(distinct, v)) 
					return (AlignedReadsSoftclip) ard.getVariation(distinct, v);
			return null;
		}
		
		private static boolean isPolya(AlignedReadsSoftclip sc, int minPrimerLen) {
			if (sc.getReadSequence().length()<minPrimerLen) return false;
			for (int i=0; i<minPrimerLen; i++)
				if (sc.getReadSequence().charAt(i)!='A') return false;
			return true;
		}

		private static ExtendedIterator<ImmutableReferenceGenomicRegion<SparseMemoryCountArray>> findPolya(
				Genomic g,
				ImmutableReferenceGenomicRegion<ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>> or, 
				int maxFrag, int minPrimerLen, double minPrimerFrac, int mergeDist, PolyaStat stat) {
			
			int numCond = or.getData().get(0).getData().getNumConditions();
			IntervalTree<UnitInterval,int[]> endCounter = new IntervalTree<>(null);
			
//			if (or.getRegion().getStart()==28466002)
//				System.out.println(or.toLocationString());
			
			// Identify all sites with >=6 A softclip
			for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read : or.getData()) {
				for (int d=0; d<read.getData().getDistinctSequences(); d++) {
					AlignedReadsSoftclip sc = getSoftclip(read.getData(), d);
					if (sc!=null && isPolya(sc,minPrimerLen)) {
						int stop = GenomicRegionPosition.ThreePrime.position(read);
						int c = read.getData().getTotalCountForDistinctInt(d, ReadCountMode.Unique);	
						endCounter.computeIfAbsent(new UnitInterval(stop), x->new int[2])[0]+=c;
					}
				}
			}

			// add the number of read-through reads
			for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read : or.getData()) {
				int c = read.getData().getTotalCountOverallInt(ReadCountMode.Unique);	
				for (int p=0; p<read.getRegion().getNumParts(); p++)
					endCounter.ei(read.getRegion().getStart(p),read.getRegion().getStop(p)).forEachRemaining(e->e.getValue()[1]+=c);
			}

			// map everything to the parent (or w/o introns)!
			ImmutableReferenceGenomicRegion<Void> parent = new ImmutableReferenceGenomicRegion<>(or.getReference(), or.getRegion().removeIntrons());
			
			IntArrayList sitesc = new IntArrayList();
			IntArrayList polyaCountsc = new IntArrayList();
			Iterator<Entry<UnitInterval, int[]>> it = endCounter.entrySet().iterator();
			while (it.hasNext()) {
				Entry<UnitInterval, int[]> en = it.next();
				if (en.getValue()[0]>=minPrimerFrac*en.getValue()[1]) {
					int p = parent.induce(en.getKey().getStart());
					sitesc.add(p);
					polyaCountsc.add(en.getValue()[0]);
				}
			}
			
			int[] sites = sitesc.toIntArray();
			int[] polyaCounts = polyaCountsc.toIntArray();
			ArrayUtils.parallelSort(sites, polyaCounts);
			SparseMemoryCountArray[] counts = new SparseMemoryCountArray[sites.length];
			for (int i=0; i<counts.length; i++) {
				counts[i] = new SparseMemoryCountArray(numCond);
			}
			
			
			IntervalTree<ArrayGenomicRegion,Transcript> hereTrans = new IntervalTree<>(null);
			
			for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(or).loop()) {
				ArrayGenomicRegion r = parent.induce(parent.getRegion().intersect(t.getRegion()));
				hereTrans.put(r,t.getData());
			}
			
			ArrayGenomicRegion non = new ArrayGenomicRegion();
			int totalNon = 0;
			
			for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read : or.getData()) {
				
				ArrayGenomicRegion rr = parent.induce(read.getRegion());
				
				int bestSite = -1;
				int bestDist = maxFrag+1;
				
				// first check unspliced
				int is = Arrays.binarySearch(sites, rr.getStop());
				if (is<0) is=-is-1;
				if (is<sites.length && sites[is]-rr.getStop()+rr.getTotalLength()<=maxFrag) {
					// found it
					bestSite = is;
					bestDist = sites[is]-rr.getStop()+rr.getTotalLength();
				} 
					
				for (ArrayGenomicRegion t : hereTrans.keys(rr.getStart(), rr.getStart()+1).filter(t->t.containsUnspliced(rr)).loop()) {
					int endont = t.induce(rr.getStop());
					int upto = Math.min(endont+maxFrag-rr.getTotalLength(),t.getTotalLength());
					ArrayGenomicRegion downstreamofread = t.map(new ArrayGenomicRegion(endont,upto));
					// this is the region where a sites is supposed to be. Find the first one!
					int iss = Arrays.binarySearch(sites, downstreamofread.getStart());
					if (iss<0) iss=-iss-1;
					int ess = Arrays.binarySearch(sites, downstreamofread.getStop());
					if (ess<0) ess=-ess-1;
					for (int i=iss; i<ess; i++) {
						if (downstreamofread.contains(sites[i])) {
							int dist = downstreamofread.induce(sites[i]);
							if (dist<bestDist) {
								bestDist = dist;
								bestSite = i;
							}
						}
					}
				}
					
				if (bestSite>=0) {
					if (stat.inpolya!=null)
						read.getData().addTotalCountsForConditions(stat.inpolya, ReadCountMode.Unique);
					read.getData().addTotalCountsForConditions(counts[bestSite], ReadCountMode.Unique);
				} else {
					if (stat.nonpolya!=null)
						read.getData().addTotalCountsForConditions(stat.nonpolya, ReadCountMode.Unique);
					non = non.union(read.getRegion());
					totalNon+=read.getData().getTotalCountOverallInt(ReadCountMode.Unique);
				}
				
			}

			if (non.getTotalLength()>0)
				stat.nonRegions.add(or.getReference(), non, new MutableInteger(totalNon));
			
			ArrayList<ImmutableReferenceGenomicRegion<SparseMemoryCountArray>> re = new ArrayList<ImmutableReferenceGenomicRegion<SparseMemoryCountArray>>();
			for (int i=0; i<sites.length; i++) {
				int e=i+1;
				int m=i;
				SparseMemoryCountArray co = counts[i];
				for (; e<sites.length && sites[e]-sites[e-1]<mergeDist; e++) {
					if (polyaCounts[e]>polyaCounts[m])
						m = e;
					co.add(counts[e]);
				}
				i = e-1;
				ArrayGenomicRegion reg = new ArrayGenomicRegion(sites[m],sites[m]+1);
				re.add(new ImmutableReferenceGenomicRegion<SparseMemoryCountArray>(or.getReference(), parent.map(reg),co));
			}
			return EI.wrap(re);
		}
	}
	
	
//	public static class PolyaClusterPositionInferenceProgram extends GediProgram {
//
//		public PolyaClusterPositionInferenceProgram(PolyaParameterSet params) {
//			addInput(params.prefix);
//			addInput(params.reads);
//			addInput(params.genomic);
//			addInput(params.nthreads);
//			addInput(params.plot);
//			
//			addOutput(params.infer3pparam);
//			addOutput(params.infer3ptsv);
//			addOutput(params.infer3pmat);
//		}
//
//		@Override
//		public String execute(GediProgramContext context) throws Exception {
//			String prefix = getParameter(0);
//			GenomicRegionStorage<DefaultAlignedReadsData> reads = getParameter(1);
//			Genomic genomic = getParameter(2);
//			int nthreads = getIntParameter(3);
//			boolean plot = getParameter(4);
//
//			
//			String[] signals = {"AATAAA","ATTAAA","AGTAAA","TATAAA"};
//			Trie<Integer> strie = new Trie<>();
//			for (int i=0; i<signals.length; i++)
//				strie.put(signals[i], i);
//			
//			int[] signalrange = {12,30}; // most of the ccds transcripts have the start of the aataaa in there! 
//			int minDistOther = 500;
//			
//			int[] count = new int[6];
//			
//			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> niceTranscripts = genomic.getTranscripts().ei()
////					EI.singleton(genomic.getTranscriptMapping().apply("ENST00000475658").toImmutable())
//				.progress(context.getProgress(), (int)genomic.getTranscripts().size(), t->t.getData().getTranscriptId())
//				.sideEffect(x->count[0]++)
////				.filter(t->t.getData().isCoding())
//				.sideEffect(x->count[1]++)
//				.filter(t->{
//					int t3p = GenomicRegionPosition.ThreePrime.position(t);
//					int t3p1 = GenomicRegionPosition.ThreePrime.position(t,1);
//					return genomic.getTranscripts()
//						.ei(new ImmutableReferenceGenomicRegion<>(t.getReference(), new ArrayGenomicRegion(Math.min(t3p1,t3p),Math.max(t3p1,t3p))))
//						.filter(o->!o.equals(t))
//						.filter(o->Math.abs(GenomicRegionPosition.ThreePrime.position(o)-t3p)<minDistOther)
//						.count()==0;
//				})
//				.filter(t->{ // same ending: take the one with the longest last exon (or the lexicographically smallest)
//					return genomic.getTranscripts()
//						.ei(t)
//						.filter(o->!o.equals(t) && GenomicRegionPosition.ThreePrime.position(o)==GenomicRegionPosition.ThreePrime.position(t))
//						.filter(o->lastExon(o).getRegion().getTotalLength()>=lastExon(t).getRegion().getTotalLength())
//						.filter(o->lastExon(o).getRegion().getTotalLength()==lastExon(t).getRegion().getTotalLength() && o.getData().getTranscriptId().compareTo(t.getData().getTranscriptId())<0)
//						.count()==0;
//				})
//				.filter(t->{ // intersects with exon intron boundary of other tr
//					return genomic.getTranscripts()
//						.ei(lastExon(t))
//						.filter(o->o.getRegion().invert().intersects(lastExon(t).getRegion()))
//						.count()==0;
//				})
//				.sideEffect(x->count[2]++)
//				.map(t->new MutableTriple<>(t,getSuffix(genomic.getSequence(lastExon(t)).toString(),minDistOther),genomic.getSequenceSave(t.getDownstream(50)).toString()))
//				.filter(tr->{
//					String expectPas = tr.Item2.substring(Math.max(0, tr.Item2.length()-signalrange[1]));
//					int paspos = expectPas.indexOf(signals[0]);
//					return paspos>=0 && paspos<signalrange[1]-signalrange[0];
//				})
//				.sideEffect(x->count[3]++)
//				.filter(tr->tr.Item2.length()>signalrange[1]?
//					strie.iterateAhoCorasick(tr.Item2.substring(0,tr.Item2.length()-signalrange[1])).count()==0
//					:true
//				)
//				.sideEffect(x->count[4]++)
//				.filter(tr->
//					SequenceUtils.getPolyAStretches(tr.Item2+tr.Item3).size()==0
//				)
//				.sideEffect(x->count[5]++)
//				.map(tr->tr.Item1)
//				.list();
//			
//			context.logf("Initial transcripts:  %d",count[0]);
//			context.logf("Removed non-coding:   %d",count[1]);
//			context.logf("Removed unclear 3':   %d",count[2]);
//			context.logf("Removed no PAS:       %d",count[3]);
//			context.logf("Removed another PAS:  %d",count[4]);
//			context.logf("Removed pA stretches: %d",count[5]);
//			
//			Collections.shuffle(niceTranscripts);
//			
//			ClusterPositionStatistics stat = EI.wrap(niceTranscripts)
//				.parallelized(nthreads, 8, ei->ei.map(tr->{
//					ImmutableReferenceGenomicRegion<Transcript> lastEx = lastExon(tr);
//					ClusterPositionStatistics tmvo = new ClusterPositionStatistics(minDistOther,tr.getData().getTranscriptId()+"\t"+(tr.getData().isCoding()?1:0));
//					for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read : reads.ei(lastEx).filter(r->lastEx.getRegion().contains(r.map(0))).loop()) {
//						int s = read.map(0);
//						if (lastEx.getRegion().contains(s)) {
//							int dist = lastEx.getRegion().getTotalLength()-lastEx.induce(s);
//							if (dist<=minDistOther) 
//								tmvo.add(dist, read.getData().getTotalCountOverall(ReadCountMode.Weight));
//						}
//					}
//					return new MutablePair<>(tmvo,tr);
//				}))
//				.progress(context.getProgress(),(int)niceTranscripts.size(), t->t.Item2.toLocationString())
//				.reduce(new ClusterPositionStatistics(),(a,b)->b.add(a.Item1));
//
//			FileUtils.writeAllLines(new String[] {"Name\tValue","Mean\t"+stat.mvo.getMean(),"Sd\t"+stat.mvo.getStandardDeviation(),""}, getOutputFile(0));
//
//			EI.seq(0, stat.histo.length).filterInt(i->stat.histo[i]>0).map(i->i+"\t"+stat.histo[i]).print("Distance\tFrequency", getOutputFile(1).getPath());
//			
//			RDataWriter rd = new RDataWriter(new FileOutputStream(getOutputFile(2)));
//			rd.writeHeader();
//			rd.write("mat", stat.getMatrix());
//			rd.write("names",stat.getNames());
//			rd.finish();
//			
//			
//			if (plot) {
//				try {
//					context.getLog().info("Running R scripts for plotting");
//					RRunner r = new RRunner(prefix+".plotclusterdist.R");
//					r.set("prefix",prefix);
//					r.addSource(getClass().getResourceAsStream("/resources/R/plotclusterdist.R"));
//					r.run(true);
//				} catch (Throwable e) {
//					context.getLog().log(Level.SEVERE, "Could not plot!", e);
//				}
//			}
//			
//			return null;
//		}
//		
//		private static String getSuffix(String s, int l) {
//			if (s.length()>l) return s.substring(s.length()-l);
//			return s;
//		}
//		
//		private static ImmutableReferenceGenomicRegion<Transcript> lastExon(ImmutableReferenceGenomicRegion<Transcript> t) {
//			GenomicRegion ex = t.getRegion().getPart(t.getReference().isMinus()?0:(t.getRegion().getNumParts()-1)).asRegion();
//			return new ImmutableReferenceGenomicRegion<>(t.getReference(), new ArrayGenomicRegion(ex),t.getData());
//		}
//	}
//
//	private static class ClusterPositionStatistics {
//		private WeightedMeanVarianceOnline mvo = new WeightedMeanVarianceOnline();
//		private double[] histo;
//		private String name;
//		
//		private ArrayList<String> names;
//		private ArrayList<double[]> all;
//		
//		public ClusterPositionStatistics(int len, String name) {
//			this.histo = new double[len+1];
//			this.name = name;
//		}
//		public double[][] getMatrix() {
//			return EI.wrap(all).toArray(double[].class);
//		}
//		
//		public String[] getNames() {
//			return EI.wrap(names).toArray(String.class);
//		}
//		public ClusterPositionStatistics() {
//			all = new ArrayList<double[]>();
//			names = new ArrayList<String>();
//		}
//		
//		public void add(int d, double w) {
//			mvo.add(d,w);
//			histo[d]+=w;
//		}
//		
//		
//		public ClusterPositionStatistics add(ClusterPositionStatistics other) {
//			this.mvo.add(other.mvo);
//			this.histo = ArrayUtils.add(this.histo, other.histo);
//			this.all.add(other.histo);
//			this.names.add(other.name);
//			return this;
//		}
//	}
	
	public static class PolyaParameterSet extends GediParameterSet {

		public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
		public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "Output prefix", true, new StringParameterType());
		public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"g", "Genomic name", true, new GenomicParameterType());
		public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "The mapped reads from the ribo-seq experiment.", false, new StorageParameterType<AlignedReadsData>());
		public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Produce plots", false, new BooleanParameterType());

		public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "Only search against chromosome 22", false, new BooleanParameterType());
		
		public GediParameter<Integer> minReads = new GediParameter<Integer>(this,"minReads", "minimal number of reads for filtered per condition", false, new IntParameterType(), 10);
		public GediParameter<Integer> minConditions = new GediParameter<Integer>(this,"minConditions", "minimal number of conditions with minReads (negative numbers are interpreted as percentage of the number of conditions)", false, new IntParameterType(), 2);
		
		
		public GediParameter<Integer> maxFrag = new GediParameter<Integer>(this,"maxFrag", "maximal fragment length to consider a read", false, new IntParameterType(), 500);
		public GediParameter<Integer> minPrimerLen = new GediParameter<Integer>(this,"minPrimerLen", "minimal A softclip length", false, new IntParameterType(), 6);
		public GediParameter<Integer> mergeDist = new GediParameter<Integer>(this,"mergeDist", "PolyA sites to merge", false, new IntParameterType(), 20);
		public GediParameter<Double> minPrimerFrac = new GediParameter<Double>(this,"minPrimerFrac", "minimal fraction of reads to infer polyA site", false, new DoubleParameterType(), 0.3);
		public GediParameter<File> readStat = new GediParameter<File>(this,"${prefix}.readstat.tsv", "Read stat file.", false, new FileParameterType());
		public GediParameter<File> nonCit = new GediParameter<File>(this,"${prefix}.outside.cit", "read clusters non associated with polya.", false, new FileParameterType());

		
		public GediParameter<File> clustercit = new GediParameter<File>(this,"${prefix}.cluster.cit", "Clusters for the viewer", false, new FileParameterType());
		public GediParameter<File> annotation = new GediParameter<File>(this,"${prefix}.annot.tsv.gz", "Annotation file", false, new FileParameterType());
		public GediParameter<File> pairwise = new GediParameter<File>(this,"${prefix}.pairwise.tsv.gz", "Pairwise file", false, new FileParameterType());
		
		public GediParameter<File> infer3ptsv = new GediParameter<File>(this,"${prefix}.3p.tsv", "Table of 3' end information.", false, new FileParameterType());
		public GediParameter<File> infer3pmat = new GediParameter<File>(this,"${prefix}.3p.Rdata", "Matrix of 3' end information.", false, new FileParameterType());
		public GediParameter<File> infer3pparam = new GediParameter<File>(this,"${prefix}.3p.parameter", "3' parameter file.", false, new FileParameterType());

	}

	
}
