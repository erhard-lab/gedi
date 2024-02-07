package gedi.slam;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.TreeMap;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;
import gedi.util.datastructure.tree.redblacktree.IntervalTreeSet;
import gedi.util.datastructure.tree.redblacktree.SimpleInterval;
import gedi.util.functions.EI;
import gedi.util.math.stat.counting.Counter;

public class NumiAlgorithm {

	private ImmutableReferenceGenomicRegion<String> gene;
	private int minEvidence = 2;
	
	public NumiAlgorithm(ImmutableReferenceGenomicRegion<String> gene) {
		this.gene = gene;
	}

	private IntArrayList posCollect = new IntArrayList();
	
	private IntervalTree<GenomicRegion, ArrayList<int[]>>[] itree = null;
	
	/**
	 * In the genomic coordinate system
	 * @param pos
	 */
	public void addConversion(int pos) {
		if (gene.getRegion().contains(pos)) {
			posCollect.add(pos);
		}
	}

	public void finishRead(GenomicRegion region, AlignedReadsData rd, int distinct) {
		if (itree==null) itree = new IntervalTree[rd.getNumConditions()];
		if (posCollect.isEmpty()) return;
		
//		region = gene.induce(region);
		int[] pattern = posCollect.toIntArray();
		Arrays.sort(pattern);
		
		for (int pp : pattern)
			if (!region.contains(pp)) 
				throw new RuntimeException("Cannot be!");
		for (int c : rd.getNonzeroCountIndicesForDistinct(distinct)) {
			if (itree[c]==null) itree[c] = new IntervalTree<GenomicRegion, ArrayList<int[]>>(null);
			
			ArrayList<int[]> list = itree[c].get(region);
			if (list==null) itree[c].put(region, list = new ArrayList<>());
			list.add(pattern);
		}
		
		posCollect.clear();
		
	}

	public int[] compute() {
		if (itree==null) return null;
		
		IntArrayList allpos = new IntArrayList();
		
		int[] re = new int[itree.length];
		for (int c=0; c<itree.length; c++) {
			if (itree[c]==null) continue;
//			System.out.println("Condition: "+c);
			for (ArrayList<int[]> v : itree[c].values())
				for (int[] a : v)
					for (int i : a)
						allpos.add(i);
			allpos.sort();
			allpos.unique();
			
			// to check for several locations!
//			TreeMap<Integer,IntervalTreeSet<SimpleInterval>> maxtree = new TreeMap<>();
			for (int start = 0; start<allpos.size(); start++) {
				for (int stop = start; stop<allpos.size(); stop++) {
					Counter<IntArrayList> clique = test(itree[c],allpos.getInt(start),allpos.getInt(stop));
					clique.elements().removeIf(k->clique.get(k, 0)<minEvidence);
					
					re[c] = Math.max(re[c], clique.elements().size());
//					maxtree.computeIfAbsent(clique.size(), x->new IntervalTreeSet<>(null)).add(new SimpleInterval(allpos.getInt(start),allpos.getInt(stop)));
//					if (!clique.elements().isEmpty()) {
//						System.out.println(gene.getReference()+":"+new ArrayGenomicRegion(allpos.getInt(start),allpos.getInt(stop)+1)+": "+clique.toSingleLine());
//					}
				}
			}
			
//			re[c] = 0;
//			// find maximal
//			while (!maxtree.isEmpty()) {
//				int tr = maxtree.lastKey();
//				IntervalTreeSet<SimpleInterval> mtr = maxtree.remove(tr);
//				int count = 0;
//				while (!mtr.isEmpty()) {
//					SimpleInterval si = mtr.first();
//					ArrayList<SimpleInterval> inter = mtr.getIntervalsIntersecting(si.getStart(), si.getEnd(), new ArrayList<>());
//					mtr.removeAll(inter);
//					count++;
//				}
////				System.out.println("Clique: "+tr+" non-overlapping count: "+count);
//				if (count>=minLocations) {
//					re[c] = tr;
//					break;
//				}
//			}
			
			allpos.clear();
		}
		
		
		return re;
	}

	private Counter<IntArrayList> test(IntervalTree<GenomicRegion, ArrayList<int[]>> itree, int start, int stop) {
//		String s = EI.wrap(itree.iterateIntervalsIntersecting(start, stop, r->r.contains(start) & r.contains(stop))).map(e->e.getKey()+"\t"+EI.wrap(e.getValue()).concat(",")).concat("\n");
		Counter<IntArrayList> set = EI.wrap(itree.iterateIntervalsIntersecting(start, stop, r->r.contains(start) & r.contains(stop)))
								.unfold(e->EI.wrap(e.getValue()))
								.map(a->restrict(a,start,stop))
								.removeNulls()
								.add(new Counter<IntArrayList>());
		return set;
	}

	private IntArrayList restrict(int[] a, int start, int stop) {
		IntArrayList re = new IntArrayList();
		for (int i=0; i<a.length; i++) {
			if (a[i]>=start && a[i]<=stop)
				re.add(a[i]);
		}
		if (re.size()==0) return null;
		return re;
	}

}
