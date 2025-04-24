package gedi.core.data.reads.subreads;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Logger;

import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsDataFactory.VarIndel;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.SingleUmiAlignedReadsData;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegion.GenomicRegionArithmetic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.BiIntConsumer;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IntToBooleanFunction;
import gedi.util.mutable.MutableInteger;
import gedi.util.sequence.DnaSequence;

/**
 * Deletions are all kept, and mismatches (and insertions) within them are discarded!
 * 
 * @author erhard
 *
 */
public class GeneSpecificUmiSenseToSubreadsConverter implements ToSubreadsConverter<BarcodedAlignedReadsData> {

	private static final int MAX_COV = 3;
	
	private String[] conditions;
	private String[] cells;
	
	private HashMap<String,Integer>[] cellIndex;
	private int cellBarcodeLength;
	
	private boolean debug = false;
	
	private String[] semantic;
	
	public GeneSpecificUmiSenseToSubreadsConverter(String[] metaDataConditions, String[] conditions, String[] cells) {
		this.conditions = conditions;
		this.cells = cells;
		
		HashMap<String, Integer> cmap = EI.wrap(metaDataConditions).indexPosition();
		if (cmap.size()!=metaDataConditions.length) throw new AssertionError("Condition names are not unique!");
		cellIndex = new HashMap[metaDataConditions.length];
		for (int i=0; i<cellIndex.length; i++)
			cellIndex[i] = new HashMap<String, Integer>();
		
		cellBarcodeLength=cells[0].length();
		if (conditions.length!=cells.length) throw new AssertionError("Conditions and cells must have equal length!");
		for (int i=0; i<conditions.length; i++) {
			Integer cind=cmap.get(conditions[i]);
			if (cind==null) throw new AssertionError("Condition name unknown: "+conditions[i]);
			if (cells[i].length()!=cellBarcodeLength) throw new AssertionError("Cell barcode length not unique!");
			cellIndex[cind].put(new DnaSequence(cells[i]).toString(), i);
		}
		
		semantic = new String[MAX_COV];
		semantic[0] = "Single read";
		semantic[1] = "Double read";
		for (int i=2; i<MAX_COV;i++)
			semantic[i]=(i+1)+"x read";
		
		
	}

	@Override
	public ToSubreadsConverter<BarcodedAlignedReadsData> setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}

	@Override
	public String[] getSemantic() {
		return semantic;
	}
	
	@Override
	public boolean isReadByRead() {
		return false;
	}

	@Override
	public void logUsedTotal(Logger logger, int used, int total) {
		ToSubreadsConverter.super.logUsedTotal(logger, used, total);
		if (total>0 && used==0) logger.severe("No UMIs used at all! Check your barcodes!");
		else if (total>0 && used<total/2) logger.warning("Very low usage ratio. Check your barcodes!");
	}
	
	
	@Override
	public ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> convert(String id,
			ReferenceSequence reference, ArrayList<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>> cluster,
			MismatchReporter reporter,BiIntConsumer usedTotal) {
		
		// TODO: Extend reporter to measure (i) number of positions in each subread cat, and (ii) specifically here, the reads per UMI statistics
		if (cluster.isEmpty()) return EI.empty();
		
		if (debug) {
			System.out.println("Convert:");
		}
		
		int numCond = cluster.get(0).getData().getNumConditions();
		
//		HashMap<GenomicRegion,MutableInteger> factories = new HashMap<>();
		int expectedReads = 0;
		
		int[][] readsPerUmi = new int[numCond][100];
		int[][] readsPerNonUmi = new int[numCond][100];
		
		ArrayList<DnaSequence> barcodes = debug?new ArrayList<>():null;
		ArrayList<GenomicRegion> regions = new ArrayList<GenomicRegion>();
		ArrayList<SubreadsAlignedReadsData> subreads = new ArrayList<SubreadsAlignedReadsData>();
		
		HashMap<DnaSequence,ArrayList<ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData>>>[] map = new HashMap[numCond];
		for (int i=0; i<numCond; i++)
			map[i] = new HashMap<DnaSequence, ArrayList<ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData>>>();

		
		// index the barcodes
		for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r : cluster) {
		
			for (int d=0; d<r.getData().getDistinctSequences(); d++) {
				if (r.getData().getMultiplicity(d)<=1) {
					if (r.getData().hasNonzeroInformation()) {
						int[] conds = r.getData().getNonzeroCountIndicesForDistinct(d);
						for (int ci=0; ci<conds.length; ci++) {
							int cond=conds[ci];
							ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> rread = new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), new SingleUmiAlignedReadsData(r.getData(), d,cond));
							for (DnaSequence dna : r.getData().getNonZeroBarcodes(d, ci)) {
								map[cond].computeIfAbsent(dna,x->new ArrayList<>()).add(rread);
							}
						}
					} else {
						for (int cond=0; cond<numCond; cond++) {
							ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> rread = new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), new SingleUmiAlignedReadsData(r.getData(), d,cond));
							for (DnaSequence dna : r.getData().getBarcodes(d, cond)) {
								map[cond].computeIfAbsent(dna,x->new ArrayList<>()).add(rread);
							}
						}
					}
				}
			}
		}
		int used = 0;
		int total = 0;
		for (int cond=0; cond<numCond; cond++) {
			for (DnaSequence barcode : map[cond].keySet()) {
				total++;
				String cb = barcode.substring(0, cellBarcodeLength);
				Integer ccond = cellIndex[cond].get(cb);
				int[][] histo = readsPerNonUmi;
				if (ccond!=null) {
					used++;
					GenomicRegion region = EI.wrap(map[cond].get(barcode)).map(r->r.getRegion()).reduce(new GenomicRegionArithmetic(),(r,a)->a.union(r)).toRegion();
	//					GenomicRegion region = EI.wrap(map.get(barcode)).map(r->r.getRegion()).reduce((r,a)->a.union(r));
	//					AlignedReadsDataFactory fac = factories.computeIfAbsent(region, x->new AlignedReadsDataFactory(cells.length).start());
//					factories.computeIfAbsent(region, x->new MutableInteger()).N++;
					SubreadsAlignedReadsData sr = createSubread(cells.length,new ImmutableReferenceGenomicRegion<>(reference, region),map[cond].get(barcode),ccond,reporter);
					
					regions.add(region);
					subreads.add(sr);
					if (debug)
						barcodes.add(barcode);
					expectedReads++;
					
					histo = readsPerUmi;
				} 
				
				histo[cond][Math.min(map[cond].get(barcode).size(),histo[cond].length-1)]++;
			}
			
			if (reporter!=null) 
				reporter.reportConvertedReadsHistogram(readsPerUmi,readsPerNonUmi);
		}
		
		ArrayList<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> re = new ArrayList<>();
		if (regions.size()>0) {
		
			// all subreads so far correspond to a single UMI!
			// sort first by region (that many reads will be generated), then by variations/subreads (that many distincts will be generated)
			// and then by the condition where the UMI is (this is nice for merging the counts)
			int[] order = ArrayUtils.order(regions.size(),(a,b)->{
				int ret = regions.get(a).compareTo(regions.get(b));
				if (ret!=0) return ret;
				ret = SubreadsAlignedReadsData.compare(subreads.get(a),0,subreads.get(b),0);
				if (ret!=0) return ret;
				return Integer.compare(subreads.get(a).getNonzeroCountIndicesForDistinct(0)[0], subreads.get(b).getNonzeroCountIndicesForDistinct(0)[0]);
			});
			
			int s = 0;
			for (int i=1; i<order.length; i++) {
				int i1 = order[i];
				int si = order[s];
				
				if (!regions.get(si).equals(regions.get(i1))) {
					re.add(new ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>(reference,regions.get(si),merge(reference,regions,subreads,barcodes,order,s,i)));
					s=i;
				}
			}
			re.add(new ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>(reference,regions.get(order[s]),merge(reference,regions,subreads,barcodes,order,s,regions.size())));
		}
		
		if (usedTotal!=null) usedTotal.accept(used, total);
		MutableInteger count = new MutableInteger();
		int uExpectedReads = expectedReads;
		return EI.wrap(re)
				.sideEffect(rere->{
					if (debug) 
						System.out.println(rere);
					count.N+=rere.getData().getTotalCountOverallInt(ReadCountMode.All);
				}).endAction(()->{
					if (debug) System.out.println();
					if (uExpectedReads!=count.N) throw new RuntimeException("Expected != Observed: "+uExpectedReads+" != "+count.N);
				});
	}

	

	/**
	 * order[s], order[s+1], ...,order[e-1] are the subreads to merge here!
	 * @param subreads
	 * @param order
	 * @param s
	 * @param e
	 * @return
	 */
	private SubreadsAlignedReadsData merge(ReferenceSequence reference,ArrayList<GenomicRegion> regions, ArrayList<SubreadsAlignedReadsData> subreads, ArrayList<DnaSequence> barcodes, int[] order, int s, int e) {
		if (e-s==1) return subreads.get(order[s]);
		
		int distincts = 1;
		for (int i=s+1; i<e; i++) {
			if (SubreadsAlignedReadsData.compare(subreads.get(order[i-1]), 0, subreads.get(order[i]), 0)<0)
				distincts++;
		}
		
		int numCond = subreads.get(order[s]).getNumConditions();;
		int[][] count = new int[distincts][];
		int[][] nonzeros = new int[distincts][];
		int[][] vars = new int[distincts][];
		CharSequence[][] indels = new CharSequence[distincts][];
		int[][] subreadGeom = new int[distincts][];
		int[][] gaps = new int[distincts][];
		
		int index = 0;
		int o = s;
		for (int i=s+1; i<=e; i++) {
			if (i==e || SubreadsAlignedReadsData.compare(subreads.get(order[i-1]), 0, subreads.get(order[i]), 0)<0) {
				
				// first take care of the umi counts!
				if (i-o==1) { // easy case, nothing to merge; just for efficiency
					int thisCond = subreads.get(order[o]).getNonzeroCountIndicesForDistinct(0)[0];
					count[index] = new int[] {1};
					nonzeros[index] = new int[] {thisCond};
				} else {
					IntArrayList tcount = new IntArrayList(3);
					IntArrayList tnonzeros = new IntArrayList(3);
					for (int ii=o; ii<i; ii++) {
						int thisCond = subreads.get(order[ii]).getNonzeroCountIndicesForDistinct(0)[0];
						if (tnonzeros.size()==0 || tnonzeros.getLastInt()!=thisCond)
							tnonzeros.add(thisCond);
						tcount.increment(tnonzeros.size()-1);
					}
					if (tcount.totalCount()>tcount.size())
						tcount.totalCount();
					count[index] = tcount.toIntArray();
					nonzeros[index] = tnonzeros.toIntArray();
				}
				
				SubreadsAlignedReadsData.setVarsAndSubreads(subreads.get(order[o]),0,vars,indels,subreadGeom,gaps,index);
				
				index++;
				o = i;
			}
		}
		
		SubreadsAlignedReadsData re = SubreadsAlignedReadsData.create(numCond, count, nonzeros, vars, indels, subreadGeom, gaps);
		
		if (debug) 
			System.out.println(" Merge:\n "+EI.seq(s,e).map(i->reference+":"+regions.get(order[i])+" "+subreads.get(order[i])+"\t"+barcodes.get(order[i])).concat("\n ")+"\n"+re);
		
		return re;
	}

	private SubreadsAlignedReadsData createSubread(int numCond, ImmutableReferenceGenomicRegion<Void> region, 
			ArrayList<ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData>> reads, int condition,
			MismatchReporter reporter) {
		
		if (debug) {
			System.out.println(region);
			EI.wrap(reads).map(r->" "+r.toString()+" "+EI.wrap(((BarcodedAlignedReadsData)r.getData().getParent()).getBarcodes(0, 0)).concat(",")).print();
			System.out.println(StringUtils.display(EI.wrap(reads).map(r->region.induce(r.getRegion()))));
		}
		
		if (reads.size()==1) {
			ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> r = reads.get(0);
			SubreadsAlignedReadsData re = SubreadsAlignedReadsData.createUnitCount(numCond, condition, new int[] {0}, new int[0], 
					r.getData().getVarIndels(0).map(vari->{
						vari = vari.reposition(region.induce(r.map(vari.getPosition())));
						if (!r.getReference().getStrand().equals(region.getReference().getStrand())) 
							vari=vari.complement();
						return vari;
					}).list()
					);
			
			if (reporter!=null) {
				for (ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> read : reads) {
					reporter.startDistinct(read, 0,(i)->false);
					for (int v=0; v<read.getData().getVariationCount(0); v++) {
						VarIndel vari = read.getData().getVarIndel(0, v);
						if (vari.isMismatch()) {
							vari=vari.reposition(region.induce(read.map(vari.getPosition())));
							if (!read.getReference().getStrand().equals(region.getReference().getStrand())) 
								vari=vari.complement();
							reporter.reportMismatch(v, false, true);
						}
					}
				}
			}
			
			if (debug) {
				System.out.println(re+"\n");
			}
			return re;
		}
		
		
		// all region part boundaries are in principle candidates for a change in the subread
		IntArrayList changeCandidates = new IntArrayList();
		HashMap<Integer,MutableInteger> covChange = new HashMap<>();
		for (ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> read : reads) {
			ArrayGenomicRegion ind = region.induce(read.getRegion());
			for (int p=0; p<ind.getNumParts(); p++) {
				changeCandidates.add(ind.getStart(p));
				covChange.computeIfAbsent(ind.getStart(p), x->new MutableInteger()).N++;
				changeCandidates.add(ind.getEnd(p));
				covChange.computeIfAbsent(ind.getEnd(p), x->new MutableInteger()).N--;
			}
		}
		
		changeCandidates.sort();
		changeCandidates.unique();
		
		if (changeCandidates.getInt(0)!=0) throw new AssertionError("Cannot be!");
		
		
		// retain only the ones where the coverage changes
		IntArrayList covPos = new IntArrayList();
		IntArrayList sub = new IntArrayList();
		IntArrayList gap = new IntArrayList();
		
		int cumu = 0;
		for (int i=0; i<changeCandidates.size(); i++) {
			int cand = changeCandidates.getInt(i);
			cumu+=covChange.get(cand).N;
			
			if (isFirstInExon(cand, region)) 
				gap.add(cand);
			
			if (i<changeCandidates.size()-1 && covChange.get(cand).N!=0) {
				int covcla = Math.min(MAX_COV, cumu)-1;
				if (sub.isEmpty() || sub.getLastInt()!=covcla) {
					// if the coverage class changes:
					if (i>0)
						sub.add(cand);
					sub.add(covcla);
				}
				covPos.add(cand);
			}
			covChange.get(cand).N=cumu;
		}
		if (cumu!=0) throw new AssertionError("Cannot be!");
		// covChange is now the coverage!
		
		// now throw together all variants from all reads, and do majority vote for each position (using the coverage to get the # reads without variant there)
		GenomicRegion allDeletion = null;
		ArrayList<VarIndel> vars = new ArrayList<>();
		for (ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> read : reads) {
			for (int v=0; v<read.getData().getVariationCount(0); v++) {
				VarIndel vari = read.getData().getVarIndel(0, v);
				if (vari.isSoftclip() || vari.isInsertion()) continue;
				if (vari.isDeletion()) {
					vari=vari.reposition(region.induce(read.map(vari.getPosition())));
					if (vari!=null) {
						GenomicRegion reg = new ArrayGenomicRegion(vari.getPosition(),vari.getPosition()+vari.getReferenceSequence().length());
						allDeletion = allDeletion==null?reg:allDeletion.union(reg);
					}
				}
				else {
					// fix problems with infrequent super long reads
					vari=vari.reposition(region.induce(read.map(vari.getPosition())));
					if (!read.getReference().getStrand().equals(region.getReference().getStrand())) 
						vari=vari.complement();
					if (vari!=null)
						vars.add(vari);
				}
			}
		}
		
		if (allDeletion!=null) 
			for (int p=0; p<allDeletion.getNumParts(); p++)
				vars.add(AlignedReadsDataFactory.createDeletion(allDeletion.getStart(p), StringUtils.repeatSequence('N', allDeletion.getLength(p)), false));
		
		Comparator<VarIndel> posComp = (a,b)->Integer.compare(a.getPosition(), b.getPosition());
		vars.sort(posComp);
		
		GenomicRegion uAllDeletion = allDeletion;
		vars = EI.wrap(vars)
			.iff(allDeletion!=null, ei->ei.filter(v->!uAllDeletion.contains(v.getPosition())))
			.multiplex(posComp, 
					list->majorityVote(list,getCoverage(list.get(0).getPosition(),covPos,covChange))
					)
			.removeNulls()
			.list();
		
		if (reporter!=null) {
			IntToBooleanFunction overFun = p->{
				return getCoverage(p,covPos,covChange)>1;
			};
			for (ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> read : reads) {
				reporter.startDistinct(read, 0,overFun);
				for (int v=0; v<read.getData().getVariationCount(0); v++) {
					VarIndel vari = read.getData().getVarIndel(0, v);
					if (vari.isMismatch()) {
						vari=vari.reposition(region.induce(read.map(vari.getPosition())));
						if (!read.getReference().getStrand().equals(region.getReference().getStrand())) 
							vari=vari.complement();
						if (vari!=null) {
							int idx = Collections.binarySearch(vars, vari,posComp);
							reporter.reportMismatch(v, overFun.applyAsInt(vari.getPosition()), idx>=0);
						}
					}
				}
			}
		}
		
		SubreadsAlignedReadsData re = SubreadsAlignedReadsData.createUnitCount(numCond, condition, sub.toIntArray(), gap.toIntArray(), vars);
		
		if (debug) {
			ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> sread = new ImmutableReferenceGenomicRegion<>(region.getReference(), region.getRegion(), re);
			System.out.println(sread+"\n");
		}

		return re;
	}


	private VarIndel majorityVote(List<VarIndel> list, int coverage) {
		if (list.size()>coverage) 
			throw new AssertionError("Cannot be!");
		
		if (list.size()*2<=coverage) return null;
		if (allSame(list)) return list.get(0);
		
		HashMap<VarIndel,MutableInteger> counter = new HashMap<>();
		for (VarIndel v : list)
			counter.computeIfAbsent(v, x->new MutableInteger()).N++;
		
		for (VarIndel v : counter.keySet()) 
			if (counter.get(v).N*2>coverage)
				return v;
		
		return null;
	}

	private boolean allSame(List<VarIndel> list) {
		for (int i=1; i<list.size(); i++)
			if (!list.get(0).equals(list.get(i)))
				return false;
		return true;
	}

	private int getCoverage(int position, IntArrayList covPos, HashMap<Integer,MutableInteger> covChange) {
		int idx = covPos.binarySearch(position);
		if (idx<0) idx=-idx-2;
		return covChange.get(covPos.getInt(idx)).N;
	}

	private boolean isFirstInExon(int cand, ImmutableReferenceGenomicRegion<Void> region) {
		return cand!=0 && cand!=region.getRegion().getTotalLength() && Math.abs(region.map(cand)-region.map(cand-1))!=1;
	}

//	private static class Merger {
//
//		private GenomicRegionArithmetic regionArithmetic;
//		private ImmutableReferenceGenomicRegion<Void> reference;
//		
//		HashMap<Integer,Counter<AlignedReadsVariation>> varsPerPosInRgr;
//		
//		public Merger() {
//		}
//		
//		public void addRegion(GenomicRegion region) {
//			regionArithmetic.union(region);
//		}
//		public void finishRegion(ReferenceSequence reference) {
//			this.reference = new ImmutableReferenceGenomicRegion<Void>(reference,regionArithmetic.toRegion());
//			regionArithmetic=null;
//		}
//		
//		public void addSubreadAndVariation(ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r) {
//			// TODO Auto-generated method stub
//			
//		}
//
//		
//
//		
//		public void addVariations(AlignedReadsDataFactory fac) {
//			for (Integer pos : varsPerPosInRgr.keySet()) {
//				// majority vote
//				Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
//				AlignedReadsVariation maxVar = ctr.getMaxElement(0);
//				int total = ctr.total()[0];
//				if (maxVar!=null && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
//					fac.addVariation(maxVar.reposition(pos));
//				}
//			}
//		}
//
//		public void add(ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r, int d) {
//			cov.add(r.getRegion(), unit);
//			int vc = r.getData().getVariationCount(d);
//			for (int v=0; v<vc; v++) {
//				AlignedReadsVariation vari = r.getData().getVariation(d, v);
//				if (vari.isDeletion() || vari.isInsertion() || vari.isMismatch()) {
//					int pos = r.map(vari.getPosition());
//					if (cov.getParentRegion().getRegion().contains(pos)) { // always except when longer than DefaultAlignedReadsData.MAX_POSITION
//						pos = cov.getParentRegion().induce(pos);
//						varsPerPosInRgr.computeIfAbsent(pos, x->new Counter<>()).add(vari);
//					}
//				}
//			}
//		}
//		
//	}
//
//	
	
	public static void main(String[] args) {
		AlignedReadsDataFactory fac1 = new AlignedReadsDataFactory(1).start();
		fac1.newDistinctSequence();
		fac1.setMultiplicity(1);
		fac1.setCount(0, 1,new DnaSequence[] {new DnaSequence("AACC")});
		ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r1 = ImmutableReferenceGenomicRegion.parse("1-:100-110",fac1.createBarcode());
		ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r2 = ImmutableReferenceGenomicRegion.parse("1-:120-130",fac1.createBarcode());
		fac1.addMismatch(4, 'A', 'G', false);
		ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> r3 = ImmutableReferenceGenomicRegion.parse("1-:100-120",fac1.createBarcode());

//		AlignedReadsDataFactory fac2 = new AlignedReadsDataFactory(1).start();
//		fac2.newDistinctSequence();
//		fac2.setMultiplicity(1);
//		fac2.setCount(0, 1);
//		fac2.setGeometry(10, 0, 10);
//		ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> rp = ImmutableReferenceGenomicRegion.parse("1:100-105|110-115|120-130",fac2.create());
		
		new GeneSpecificUmiSenseToSubreadsConverter(new String[] {"C"}, new String[] {"C"}, new String[] {"AA"})
			.setDebug(true)
			.convert("Test",Chromosome.obtain("1+"), new ArrayList<>(Arrays.asList(r1,r2,r3)), null,null)
			.drain();
		
		ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> s1 = ImmutableReferenceGenomicRegion.parse("1-:100-110",new SingleUmiAlignedReadsData(r1.getData(),0,0));
		ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> s2 = ImmutableReferenceGenomicRegion.parse("1-:120-130",new SingleUmiAlignedReadsData(r2.getData(),0,0));
		ImmutableReferenceGenomicRegion<SingleUmiAlignedReadsData> s3 = ImmutableReferenceGenomicRegion.parse("1-:100-120",new SingleUmiAlignedReadsData(r3.getData(),0,0));

		SubreadsAlignedReadsData ss = new GeneSpecificUmiSenseToSubreadsConverter(new String[] {"C"}, new String[] {"C"}, new String[] {"AA"}).
				createSubread(0, ImmutableReferenceGenomicRegion.parse("1+:100-120"),new ArrayList<>(Arrays.asList(s3)), 0,null);
		
		System.out.println(ss);
		
//		System.out.println();
		
//		System.out.println(new PairedEndToSubreadsConverter(true).setDebug(true).convert(rp, true, null));
		
		
	}
	
	
}
