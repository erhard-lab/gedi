package gedi.core.data.reads;

import java.util.Arrays;
import java.util.Iterator;

import org.h2.util.IntIntHashMap;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.FunctorUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

public class AlignedReadsDataMerger {
	
	private int[] offsets;
	private AlignedReadsDataFactory fac;
	
	private IntIntHashMap[] idMapping;
	private int nextId = 1; // start with 1 as IntIntHashMap does not support 0->0
	
//	public AlignedReadsDataMerger(){
//	}
	
	public AlignedReadsDataMerger(int... numConditions) {
		fac = new AlignedReadsDataFactory(ArrayUtils.sum(numConditions));
		offsets = new int[numConditions.length+1];
		for (int i=1; i<offsets.length; i++) 
			offsets[i] = offsets[i-1]+numConditions[i-1];
		idMapping = new IntIntHashMap[numConditions.length];
		for (int i=0; i<idMapping.length; i++)
			idMapping[i] = new IntIntHashMap();
		
	}


	public ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> merge(ImmutableReferenceGenomicRegion<? extends AlignedReadsData>... data) {
		AlignedReadsData[] d2 = new AlignedReadsData[data.length];
		ReferenceSequence ref = null;
		GenomicRegion region = null;
		for (int i=0; i<d2.length; i++) {
			if (data[i]!=null) {
				d2[i] = data[i].getData();
				if (ref==null) ref = data[i].getReference();
				else if (!ref.equals(data[i].getReference())) throw new IllegalArgumentException("ReferenceRegions must be equal!");
				if (region==null) region = data[i].getRegion();
				else if (!region.equals(data[i].getRegion())) throw new IllegalArgumentException("ReferenceRegions must be equal!");
			}
		}
		return new ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>(ref,region,merge(d2));
	}
	
	private static class DistinctSeq implements Comparable<DistinctSeq> {
		int d;
		AlignedReadsVariation[] var;
		int geom;
		public DistinctSeq(int d, AlignedReadsVariation[] var, int geom) {
			this.d = d;
			this.var = var;
			this.geom = geom;
			Arrays.sort(var);
		}
		@Override
		public int compareTo(DistinctSeq o) {
			int n = Math.min(var.length,o.var.length);
			for (int i=0; i<n; i++) {
				int r = var[i].compareTo(o.var[i]);
				if (r!=0)
					return r;
			}
			int re = var.length-o.var.length;
			if (re==0) {
				return Integer.compare(geom, o.geom);
			}
			return re;
		}
	}

	public DefaultAlignedReadsData merge(AlignedReadsData... data) {
//		if (fac==null) {
//			fac = new AlignedReadsDataFactory(data.length);
//			offsets = new int[data.length+1];
//			for (int i=1; i<offsets.length; i++)
//				offsets[i] = i;
//		}
		fac.start();
		
		int max = 1;
		for (int i=0; i<data.length; i++)
			if (data[i]!=null)
				max = Math.max(max,data[i].getDistinctSequences());
		
		
		Iterator<DistinctSeq>[] it = new Iterator[data.length];
		for (int i=0; i<data.length; i++) 
			if (data[i]!=null) {
				DistinctSeq[] vars = new DistinctSeq[data[i].getDistinctSequences()];
				for (int v=0; v<data[i].getDistinctSequences(); v++) 
					vars[v] = new DistinctSeq(v,data[i].getVariations(v),data[i].getRawGeometry(v));
				Arrays.sort(vars);
				it[i] = EI.wrap(vars);
			} else
				it[i] = EI.empty();
		
		
		ExtendedIterator<DistinctSeq[]> pit = EI.parallel(DistinctSeq.class, FunctorUtils.naturalComparator(), it);
		while (pit.hasNext()) {
			DistinctSeq[] n = pit.next();
			
			int i;
			for (i=0;n[i]==null; i++);
			int v = n[i].d;
			
			fac.newDistinctSequence();
			if (data[i].hasId())
				fac.setId(getId(i,data[i].getId(v), data[i].getMultiplicity(v)>1));
			fac.setMultiplicity(data[i].getMultiplicity(v));
			if (data[i].hasWeights())
				fac.setWeight(data[i].getWeight(v));
			if (data[i].hasGeometry())
				fac.setGeometry(data[i].getGeometryBeforeOverlap(v), data[i].getGeometryOverlap(v),data[i].getGeometryAfterOverlap(v));
			for (AlignedReadsVariation vari : n[i].var)
				fac.addVariation(vari);
			
			BarcodedAlignedReadsData bc = data[i] instanceof BarcodedAlignedReadsData?(BarcodedAlignedReadsData)data[i]:null;
			
			if (data[i].hasNonzeroInformation()) {
				int[] cs = data[i].getNonzeroCountIndicesForDistinct(v);
				for (int ci=0; ci<cs.length; ci++) {
					int c=cs[ci];
					fac.setCount(offsets[i]+c, data[i].getNonzeroCountValueForDistinct(v, ci),bc==null?null:bc.getNonZeroBarcodes(v, ci));
				}
			}
			else {
				for (int c=0; c<data[i].getNumConditions(); c++)
					if (data[i].getCount(v, c)>0)
						fac.setCount(offsets[i]+c, data[i].getCount(v, c), bc==null?null:bc.getBarcodes(v, c));
			}
			
			for (i++ ; i<n.length; i++) {
				bc = data[i] instanceof BarcodedAlignedReadsData?(BarcodedAlignedReadsData)data[i]:null;
				if (n[i]!=null) {
					if (data[i].hasNonzeroInformation()) {
						int[] cs = data[i].getNonzeroCountIndicesForDistinct(n[i].d);
						for (int ci=0; ci<cs.length; ci++) {
							int c=cs[ci];
							fac.setCount(offsets[i]+c, data[i].getNonzeroCountValueForDistinct(n[i].d, ci),bc==null?null:bc.getNonZeroBarcodes(n[i].d, ci));
						}
					} 
					else {
						for (int c=0; c<data[i].getNumConditions(); c++)
							if (data[i].getCount(n[i].d, c)>0)
								fac.setCount(offsets[i]+c, data[i].getCount(n[i].d, c),  bc==null?null:bc.getBarcodes(n[i].d, c));	
					}
				}
			}
			
//			if (barcodes) {
//				for (i=0; i<n.length; i++) {
//					if (n[i]!=null) 
//						for (int c=0; c<data[i].getNumConditions(); c++)
//							if (((BarcodedAlignedReadsData)data[i]).getBarcodes(n[i].d, c).length>0)
//								fac.addBarcode(offsets[i]+c, ((BarcodedAlignedReadsData)data[i]).getBarcodes(n[i].d, c));
//				}
//			}
			
		}
//		
//		BitVector done = new BitVector(data.length*max);
//		// do a linear search, as only few variations are expected! HAHAHA, you fool! is 1231723520 "few"???
//		for (int i=0; i<data.length; i++)
//			if (data[i]!=null) {
//				for (int v=0; v<data[i].getDistinctSequences(); v++) {
//					if (done.getQuick(i+v*data.length)) continue;
//					fac.newDistinctSequence();
//					fac.setId(getId(i,data[i].getId(v), data[i].getMultiplicity(v)>1));
//					fac.setMultiplicity(data[i].getMultiplicity(v));
//					if (data[i].hasWeights())
//						fac.setWeight(data[i].getWeight(v));
//					if (data[i].hasGeometry())
//						fac.setGeometry(data[i].getGeometryBeforeOverlap(v), data[i].getGeometryOverlap(v),data[i].getGeometryAfterOverlap(v));
//					for (int c=0; c<data[i].getNumConditions(); c++)
//						fac.setCount(offsets[i]+c, data[i].getCount(v, c));
//					
//					AlignedReadsVariation[] vars = data[i].getVariations(v);
//					Arrays.sort(vars);
//					
//					for (AlignedReadsVariation vari : vars)
//						fac.addVariation(vari);
//					
//					for (int j=i+1; j<data.length; j++) 
//						if (data[j]!=null) {
//							for (int w=0; w<data[j].getDistinctSequences(); w++) {
//								if (done.getQuick(j+w*data.length)) continue;
//								
//								AlignedReadsVariation[] wars = data[j].getVariations(w);
//								Arrays.sort(wars);
//								
//								if (Arrays.equals(vars, wars)) {
//									for (int c=0; c<data[j].getNumConditions(); c++)
//										fac.setCount(offsets[j]+c, data[j].getCount(w, c));
//									done.putQuick(j+w*data.length,true);
//									break;
//								}
//							}
//						}
//					
//				}
//			}
		
		DefaultAlignedReadsData re = fac.createDefaultOrBarcode();
//		System.out.println(Arrays.toString(data)+" -> "+re);
		return re;
	}


	private int getId(int file, int oldId, boolean save) {
		int id = idMapping[file].get(oldId);
		if (id==IntIntHashMap.NOT_FOUND) {
			if (save)
				idMapping[file].put(oldId, nextId);
			id = nextId++;
		}
		return id;
	}
	
}
