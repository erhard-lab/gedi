package gedi.grand10x;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand10x.javapipeline.Grand10XDemultiplexProgram.SnpData;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.genomic.CoverageAlgorithm;
import gedi.util.math.stat.counting.Counter;

public class Grand10XMerger {

		
		private static NumericArray unit = NumericArray.wrap(1);
		HashMap<Integer,Counter<AlignedReadsVariation>> varsPerPosInRgr;
		CoverageAlgorithm cov;
		
		GenomicRegion region;
		
		public Grand10XMerger() {
		}
		
		public ReadCategory classify(Genomic g, TrimmedGenomicRegion buf) {
			return ReadCategory.classify(g, cov.getParentRegion().getReference(), buf.set(cov.getParentRegion().getRegion()));
		}
		
		public void computeStatistics(MismatchPerCoverageStatistics stat, CharSequence sequence, String cat, TreeMap<Integer, SnpData> snpdata) {
				
			for (int i=0; i<cov.getParentRegion().getRegion().getTotalLength(); i++)  
				if (!snpdata.containsKey(cov.getParentRegion().map(i)))
					stat.addCoverage(cat,sequence.charAt(i), cov.getCoverages(i).getInt(0));
				
			for (Integer pos : varsPerPosInRgr.keySet()) {
				if (snpdata.containsKey(cov.getParentRegion().map(pos)))
					continue;
				// majority vote
				Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
				AlignedReadsVariation maxVar = ctr.getMaxElement(0);
				int total = ctr.total()[0];
				if (maxVar!=null && maxVar.isMismatch() && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
					stat.addMismatch(cat,maxVar.getReferenceSequence().charAt(0), maxVar.getReadSequence().charAt(0), cov.getCoverages(pos).getInt(0));
				}
			}
		}
		
		
		public void addRegion(GenomicRegion reg) {
			if (region==null) region = reg;
			else region = region.union(reg);
		}
		
		public void finishedRegion(ReferenceSequence ref) {
			cov = new CoverageAlgorithm(ref,region);
			varsPerPosInRgr = new HashMap<>();
		}
		
		public GenomicRegion getRegionWithCoverageMin(int mincov) {
			NumericArray prof = cov.getProfile(0);
			int start = -1; 
			int stop = -1;
			for (int i=0; i<prof.length(); i++) {
				if (start==-1) {
					if (prof.getInt(i)>=mincov) 
						start = stop = i;
				}
				else {
					if (prof.getInt(i)>=mincov) 
						stop = i;
				}
			}
			if (start==-1) return null;
			return cov.getParentRegion().map(new ArrayGenomicRegion(start, stop+1));
		}
		
		public void addSnpData(TreeMap<Integer, SnpData> snpdata) {
			for (Integer pos : varsPerPosInRgr.keySet()) {
				// majority vote
				Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
				AlignedReadsVariation maxVar = ctr.getMaxElement(0);
				int total = ctr.total()[0];
				if (maxVar!=null && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
					snpdata.computeIfAbsent(cov.getParentRegion().map(pos), x->new SnpData()).incrementMM();
				}
			}
		}
		
		public void addSnpCoverage(TreeMap<Integer, SnpData> snpdata) {
			Iterator<Entry<Integer, SnpData>> it = snpdata.subMap(cov.getParentRegion().getRegion().getStart(), true, cov.getParentRegion().getRegion().getEnd(), false).entrySet().iterator();
			while (it.hasNext()) {
				Entry<Integer, SnpData> en = it.next();
				en.getValue().incrementCov();
			}
		}
		
		public void addVariations(GenomicRegion target, AlignedReadsDataFactory fac) {
			target = cov.getParentRegion().induce(target);
			
			for (Integer pos : varsPerPosInRgr.keySet()) {
				
				if (target.contains(pos)) {
					// majority vote
					Counter<AlignedReadsVariation> ctr = varsPerPosInRgr.get(pos);
					AlignedReadsVariation maxVar = ctr.getMaxElement(0);
					int total = ctr.total()[0];
					if (maxVar!=null && ctr.get(maxVar)[0]>cov.getCoverages(pos).getInt(0)-total) {
						fac.addVariation(maxVar.reposition(target.induce(pos)));
					}
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

		public ImmutableReferenceGenomicRegion<?> getParentRegion() {
			return cov.getParentRegion();
		}

		
		
	}