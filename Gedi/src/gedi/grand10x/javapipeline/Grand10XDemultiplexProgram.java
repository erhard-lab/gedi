package gedi.grand10x.javapipeline;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.BiFunction;
import java.util.function.Supplier;

import executables.Grand10X;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.grand10x.MismatchPerCoverageStatistics;
import gedi.grand10x.Grand10XMerger;
import gedi.grand10x.ReadCategory;
import gedi.grand10x.TrimmedGenomicRegion;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.ParallelizedIterator;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.StreamLineWriter;
import gedi.util.math.stat.counting.Counter;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.sequence.DnaSequence;
import jdistlib.Beta;

public class Grand10XDemultiplexProgram extends GediProgram {



		public Grand10XDemultiplexProgram(Grand10xParameterSet params) {
			addInput(params.prefix);
			addInput(params.reads);
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.maxdist);
			addInput(params.whitelist);
			addInput(params.test);
			addInput(params.infer3pparam);
			addInput(params.mincov);
			addInput(params.clustercit);
			addInput(params.computeStat);
			addInput(params.snpConv);
			addInput(params.snpPval);
			

			addOutput(params.umicit);
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
			int mincov = getIntParameter(8);
			File clustFile = getParameter(9);
			boolean stat = getBooleanParameter(10);
			double snpconv = getDoubleParameter(11);
			double snpp = getDoubleParameter(12);
			
			
			String testloc="9+";
//			testloc = "8-:71596058-71596227";
			
			double[] meansd = {Double.NaN, Double.NaN};
			HeaderLine h = new HeaderLine();
			EI.lines(pparam).header(h).split("\t").forEachRemaining(a->{
				if (a[0].equals("Mean"))
					meansd[0] = Double.parseDouble(a[1]);
				else if (a[0].equals("Sd"))
					meansd[1] = Double.parseDouble(a[1]);
			});
			
			context.logf("Mean=%.1f, Sd=%.1f", meansd[0],meansd[1]);
			
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
			
			MemoryIntervalTreeStorage<Transcript> transcripts = genomic.getTranscripts();
			
			BiFunction<ReferenceGenomicRegion<?>,TrimmedGenomicRegion,ExtendedIterator<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> readSupp = (rgr,reuse)->
				reads.ei(rgr).filter(r->{
					if (r.getData().getMultiplicity(0)>1) 
						return false;
					if (!r.getRegion().intersects(rgr.getRegion()))
						return false;
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> ctrans = transcripts.ei(r)
								.chain(transcripts.ei(r).map(t->t.toMutable().toOppositeStrand().toImmutable()))
								.filter(t->reuse.set(r.getRegion()).isCompatibleWith(t.getRegion())).list();
					
					if (r.getRegion().getNumParts()>1 && ctrans.isEmpty())
						return false;
					
					if (!rgr.getRegion().contains(r.getRegion())) {
						System.out.println();
						System.out.println("Region: "+rgr);
						System.out.println(r);
						throw new RuntimeException();
					}
					
					return true;
				});
				
				
			String[] cats = EI.wrap(ReadCategory.values())
					.unfold(rc->EI.wrap(genomic.getOrigins()).chain(EI.singleton("Mito")).map(g->g+"\t"+rc.toString()))
					.unfold(rc->EI.wrap(new String[] {"T->C","A->G","Weak"}).map(g->rc.toString()+"\t"+g))
					.toArray(String.class);
			int readlen = reads.ei().filter(r->r.getData().getVariationCount(0)==0).first().getRegion().getTotalLength();
			Supplier<MismatchPerCoverageStatistics> statCompMaker = stat?()->new MismatchPerCoverageStatistics(readlen,cats,"Genome\tCategory\tPredominant Conversion"):()->null; 
				
			CenteredDiskIntervalTreeStorage<MutableInteger> clusterCit = new CenteredDiskIntervalTreeStorage<>(clustFile.getPath());
				
			ParallelizedIterator<ImmutableReferenceGenomicRegion<MutableInteger>, ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>, MismatchPerCoverageStatistics> readit = (test?clusterCit.ei(testloc):clusterCit.ei())
				.progress(context.getProgress(), (int)clusterCit.size(), x->"Demultiplexing reads "+x.toLocationString())
				.parallelized(nthreads, 1, statCompMaker,(xei,statComp)->
					xei.unfold(x->{
						//ImmutableReferenceGenomicRegion<Void>[] peaks = findPeaks(x, readSupp, conditions, meansd[0], meansd[1], genomic);
						//if (peaks.length==0) return null;
						return demultiplex(x, readSupp, conditions, genomic, whitelists, mincov,statComp, snpconv, snpp);
		//				return EI.<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>>empty();
					}
					))
				;
				
			out.fill(readit.removeNulls(),context.getProgress());
			
			MismatchPerCoverageStatistics statComp = readit.drainStates().removeNulls().reduce((a,b)->a.add(b));
			
			if (statComp!=null) {
				LineWriter sout = new LineOrientedFile(prefix+".mismatchcov.tsv").write();
				statComp.writeCov(sout);
				sout.close();
				
				sout = new LineOrientedFile(prefix+".mismatchpos.tsv").write();
				statComp.writePos(sout);
				sout.close();
				
				sout = new LineOrientedFile(prefix+".snpdata").write();
				statComp.writeSnps(sout,snpconv);
				sout.close();
				
				sout = new LineOrientedFile(prefix+".segregating.tsv").write();
				statComp.writeSeg(sout);
				sout.close();
			}
			

			return null;
		}
		
		private static ExtendedIterator<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> demultiplex(
				ImmutableReferenceGenomicRegion<?> or,
				BiFunction<ReferenceGenomicRegion<?>,TrimmedGenomicRegion,ExtendedIterator<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>>> reads,
				String[] conditions, Genomic g, HashMap<String,Integer>[] whitelists, int mincov,
				MismatchPerCoverageStatistics stat, double snpconv, double snpp) {
			
			String locs = or.toLocationString();
			TrimmedGenomicRegion reuse2 = new TrimmedGenomicRegion();
			int cellBarcodeLength = whitelists[0].keySet().iterator().next().length();
			int nind = EI.wrap(whitelists).unfold(m->EI.wrap(m.values())).mapToInt(x->x).max()+1;
			
			
			MismatchPerCoverageStatistics checker = new MismatchPerCoverageStatistics(stat.getReadLen(), new String[] {locs}, "Cluster");
			HashMap<GenomicRegion,ArrayList<DefaultAlignedReadsData>> pre = new HashMap<>();
			
			int umis = 0;
			
			for (int c=0; c<conditions.length; c++) {
				// determine region for each bc
				HashMap<DnaSequence,Grand10XMerger> regions = new HashMap<>();
				for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : reads.apply(or,reuse2).loop()) {
					for (int d=0; d<r.getData().getDistinctSequences(); d++) {
						for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
							String cellbc = bc.toString().substring(0,cellBarcodeLength);
							Integer cell = whitelists[c].get(cellbc); 
							if (cell!=null) {
								Grand10XMerger p = regions.computeIfAbsent(bc, x->new Grand10XMerger());
								p.addRegion(r.getRegion());
							}
						}
					}
				}
				for (Grand10XMerger m : regions.values())
					m.finishedRegion(or.getReference());
				umis+=regions.size();
				
				// add to merger
				for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : reads.apply(or,reuse2).loop()) {
					for (int d=0; d<r.getData().getDistinctSequences(); d++) {
						for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
							String cellbc = bc.toString().substring(0,cellBarcodeLength);
							Integer cell = whitelists[c].get(cellbc); 
							if (cell!=null) {
								Grand10XMerger p = regions.get(bc);
								p.add(r,d);
							}
						}
					}
				}
				
				// estimate SNPs (for stats, they are added to the umiread anyway!)
				TreeMap<Integer,SnpData> snpdata = new TreeMap<Integer,SnpData>();
				for (DnaSequence pp : regions.keySet()) 
					regions.get(pp).addSnpData(snpdata);
				for (DnaSequence pp : regions.keySet()) 
					regions.get(pp).addSnpCoverage(snpdata);
				
				Iterator<SnpData> it = snpdata.values().iterator();
				while (it.hasNext()) {
					SnpData s = it.next();
					if (s.getPval(snpconv)>snpp)
						it.remove();
				}
				stat.addSnps(or.getReference().getName(),snpdata);
				
				
				// merge
				for (DnaSequence pp : regions.keySet()) {
					
					Grand10XMerger p = regions.get(pp);
					p.computeStatistics(checker, g.getSequence(p.getParentRegion()),locs,snpdata);
					
					String cellbc = pp.substring(0,cellBarcodeLength);
					Integer cell = whitelists[c].get(cellbc); 
					
					GenomicRegion reg = p.getRegionWithCoverageMin(mincov);
					if (reg!=null) {
						AlignedReadsDataFactory fac = new AlignedReadsDataFactory(nind);
						fac.start();
						fac.newDistinctSequence();
						fac.setMultiplicity(1);
						fac.setCount(cell, 1);
						p.addVariations(reg,fac);
						DefaultAlignedReadsData ard = fac.create();
						
						ArrayList<DefaultAlignedReadsData> l = pre.computeIfAbsent(reg, x->new ArrayList<>());
						l.add(ard);
						
						if (l.size()>10) {
							fac = new AlignedReadsDataFactory(nind);
							fac.start();
							for (DefaultAlignedReadsData lard : l) 
								fac.add(lard);
							fac.makeDistinct();
							l.clear();
							l.add(fac.create());
						}
					}
					
				}
				
				String pred = checker.getMmFraction(locs, 'T', 'C')<checker.getMmFraction(locs, 'A', 'G')?"A->G":"T->C";
				if (umis<100) pred="Weak";
				
				if (stat!=null) 
					for (DnaSequence pp : regions.keySet()) {
						Grand10XMerger p = regions.get(pp);
						String ori = g.getOrigin(p.getParentRegion().getReference()).toString();
						if (p.getParentRegion().getReference().isMitochondrial()) ori = "Mito";
						p.computeStatistics(stat, g.getSequence(p.getParentRegion()),ori+"\t"+p.classify(g,reuse2)+"\t"+pred,snpdata);
					}
				
				
				// per read statistics
				if (stat!=null) {
					for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : reads.apply(or,reuse2).loop()) {
						String ori = g.getOrigin(r.getReference()).toString();
						if (r.getReference().isMitochondrial()) ori = "Mito";
						
						for (int d=0; d<r.getData().getDistinctSequences(); d++) {
							for (DnaSequence bc : r.getData().getBarcodes(d, c)) {
								String cellbc = bc.toString().substring(0,cellBarcodeLength);
								Integer cell = whitelists[c].get(cellbc); 
								if (cell!=null) {
									Grand10XMerger p = regions.get(bc);
									ReadCategory cat = p.classify(g,reuse2);
									CharSequence seq = g.getSequence(r);
									int vc = r.getData().getVariationCount(d);
									for (int v=0; v<vc; v++) {
										AlignedReadsVariation vari = r.getData().getVariation(d, v);
										if (vari.isMismatch()) {
											int pos = vari.getPosition();
											if (!snpdata.containsKey(r.map(pos)))
												stat.addReadMismatch(ori+"\t"+cat+"\t"+pred, vari.getReferenceSequence().charAt(0), vari.getReadSequence().charAt(0), pos);
										}
									}
									for (int i=0; i<r.getRegion().getTotalLength(); i++)  
										if (!snpdata.containsKey(r.map(i)))
											stat.addReadCoverage(ori+"\t"+cat+"\t"+pred,seq.charAt(i), i);
								}
								
								
							}
						}
					}
				}
				
			}
			
			
			if (stat!=null && umis>=100 && 
					(checker.getMmFraction(locs, 'T', 'C')>0.01 || checker.getMmFraction(locs, 'A', 'G')>0.01)) {
				stat.addSegregatingAnalysis(locs,ReadCategory.classify(g, or.getReference(), reuse2.set(or.getRegion())),checker.getMmFraction(locs, 'T', 'C'),checker.getMmFraction(locs, 'A', 'G'),umis);
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
								fac.add(ard);
							fac.makeDistinct();
							return new ImmutableReferenceGenomicRegion<>(or.getReference(),reg,fac.create());
						}
					});

			
		}
		
		
		
		public static class SnpData {

			private int cov;
			private int mm;
			
			public SnpData() {
			}

			public double getPval(double snpconv) {
				return Beta.cumulative(snpconv, mm+snpconv, cov-mm+1, true, false);
			}

			public int getCov() {
				return cov;
			}

			public int getMm() {
				return mm;
			}

			public void incrementMM() {
				mm++;
			}

			public void incrementCov() {
				cov++;
			}

			public SnpData add(SnpData other) {
				this.cov+=other.cov;
				this.mm+=other.mm;
				return this;
			}
			
		}

	}