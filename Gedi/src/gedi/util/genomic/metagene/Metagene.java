package gedi.util.genomic.metagene;

import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.function.Consumer;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;

import gedi.core.data.annotation.Transcript;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.IndexDoubleProcessor;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.functions.EI;
import gedi.util.functions.ParallelizedState;
import gedi.util.r.RDataWriter;

public class Metagene implements ParallelizedState<Metagene>{

	private int width;
	private MetageneRange[] ranges;
	
	private ArrayList<String> names = new ArrayList<String>();
	private ArrayList<double[]> profiles = new ArrayList<double[]>();
	private ArrayList<int[]> lengths = new ArrayList<int[]>();
	
	private Consumer<double[]> normalizer = a->{};
	
	Metagene(int width, Consumer<double[]> normalizer, ArrayList<MetageneDefinition.MetageneRangeDefinition> f) {
		this.width = width;
		this.normalizer = normalizer;
		ranges = (MetageneRange[]) Array.newInstance(MetageneRange.class,f.size());
		double twidth = EI.wrap(f).mapToDouble(d->d.width).sum();
		int lastEnd = 0;
		double cumsumwidth = 0;
		for (int i=0; i<ranges.length; i++) {
			cumsumwidth+=f.get(i).width;
			int end = (int) (width*cumsumwidth/twidth);
			ranges[i] = new MetageneRange(end-lastEnd,f.get(i));
			lastEnd = end;
		}
	}

	Metagene(int width, MetageneRange[] ranges, Consumer<double[]> normalizer) {
		this.width = width;
		this.ranges = ranges;
		this.normalizer = normalizer;
	}

	public void add(String name, ImmutableReferenceGenomicRegion<?> g) {
		double[] profile = new double[width];
		int[] lengths = new int[ranges.length];
		int start = 0;
		
		for (int i=0; i<ranges.length; i++) {
			ReferenceSequence ref = ranges[i].reg.getReference(g);
			GenomicRegion reg = ranges[i].reg.getRegion(g);
			lengths[i] = reg.getTotalLength();
			
			IndexDoubleProcessor data = ranges[i].provider.getData(new ImmutableReferenceGenomicRegion<>(ref, reg));
			
			if (ranges[i].maxLength>0) {
				// truncate
				if (data.length()>ranges[i].maxLength) {
					int startInData = (int) ((data.length()-ranges[i].maxLength)*ranges[i].where);
					scale(data, startInData, startInData+ranges[i].maxLength, profile, start, start+ranges[i].width,0,ranges[i].shrinkStat);
//					if (cov!=null) {
//						scale(cov, startInData, startInData+ranges[i].maxLength, profile, start, start+ranges[i].width,1,negate(ranges[i].shrinkStat));
//						for (int ii=start+1; ii<start+ranges[i].width; ii++)
//							profile[ii] = profile[ii-1]+profile[ii];
//					}
					profile[start+ranges[i].width-1] = Double.POSITIVE_INFINITY;
				}
				// pad
				else {
					double factor = ranges[i].width/(double)ranges[i].maxLength;
					int scaledDataLen = (int) (data.length()*factor);
					int startInRe = (int) ((ranges[i].width-scaledDataLen)*ranges[i].where);
					Arrays.fill(profile, start, start+startInRe, Double.NaN);
					Arrays.fill(profile, (int)(start+startInRe+scaledDataLen), start+ranges[i].width, Double.NaN);
					scale(data, 0, data.length(), profile, start+startInRe, (int)(start+startInRe+scaledDataLen),0,ranges[i].shrinkStat);
//					if (cov!=null) {
//						scale(cov, 0, data.length(), profile, start+startInRe, (int)(start+startInRe+scaledDataLen),1,negate(ranges[i].shrinkStat));
//						for (int ii=start+startInRe+1; ii<(int)(start+startInRe+scaledDataLen); ii++)
//							profile[ii] = profile[ii-1]+profile[ii];
//					}
					
				}
			}
			else {
				if (ranges[i].fixedFactor>0)
					checkFixedFactor(0, data.length(), start,start+ranges[i].width, ranges[i].fixedFactor);
				scale(data, 0, data.length(), profile, start, start+ranges[i].width,0,ranges[i].shrinkStat);
//				if (cov!=null) {
//					scale(cov, 0, data.length(), profile, start, start+ranges[i].width,1,negate(ranges[i].shrinkStat));
//					for (int ii=start+1; ii<start+ranges[i].width; ii++)
//						profile[ii] = profile[ii-1]+profile[ii];
//				}
				
			}
			
			start+=ranges[i].width;
		}
		
		normalizer.accept(profile);
		
		this.profiles.add(profile);
		this.lengths.add(lengths);
		this.names.add(name);
		
	}
	
//	private UnivariateStatistic negate(UnivariateStatistic s) {
//		return new UnivariateStatistic() {
//			@Override
//			public double evaluate(double[] values, int begin, int length) throws MathIllegalArgumentException {
//				return -s.evaluate(values, begin, length);
//			}
//			
//			@Override
//			public double evaluate(double[] values) throws MathIllegalArgumentException {
//				return -s.evaluate(values);
//			}
//			
//			@Override
//			public UnivariateStatistic copy() {
//				return negate(s.copy());
//			}
//		};
//	}

	public void writeRds(String path) throws IOException {
		RDataWriter out = new RDataWriter(new GZIPOutputStream(new FileOutputStream(path)));
		out.writeHeader(false);
		LinkedHashMap<String, Object> attr = new LinkedHashMap<>();
		attr.put("dimnames", Arrays.asList(names.toArray(new String[0]),null));
		attr.put("rangenames", EI.wrap(ranges).map(r->r.name).toArray(String.class));
		attr.put("rangewidth", EI.wrap(ranges).mapToInt(r->r.width).toIntArray());
		attr.put("realstart", EI.wrap(ranges).map(r->r.reg.getStart()).toArray(String.class));
		attr.put("realstop", EI.wrap(ranges).map(r->r.reg.getStop()).toArray(String.class));
		attr.put("realunit", EI.wrap(ranges).map(r->r.reg.getUnit()).toArray(String.class));
		out.write(null, profiles,attr, new String[] {"matrix","metagene",});
		out.finish();
	}
	
	public void writeTsv(String path) throws IOException {
		EI.seq(0, names.size()).map(i->names.get(i)+"\t"+StringUtils.concat("\t", profiles.get(i))).print(path);
	}
	
	public void writeLengths(String path) throws IOException {
		String header = "Item\t"+EI.wrap(ranges).map(r->r.name).concat("\t");
		EI.seq(0, names.size()).map(i->names.get(i)+"\t"+StringUtils.concat("\t", lengths.get(i))).print(header,path);
	}
	
	@Override
	public Metagene spawn(int index) {
		return new Metagene(width, ranges,normalizer);
	}


	@Override
	public void integrate(Metagene other) {
		names.addAll(other.names);
		profiles.addAll(other.profiles);
		lengths.addAll(other.lengths);
	}	

	public double[] getProfile(String name) {
		return profiles.get(names.indexOf(name));
	}
	
	private static void checkFixedFactor(int dstart, int dend, int start, int end, double fixedFactor) {
		if (fixedFactor>0 && (int)Math.round((end-start)/fixedFactor)!=dend-dstart)
			throw new RuntimeException("Not the fixed factor: "+dstart+","+dend+","+start+","+end+","+fixedFactor);
	}
	
	private static void scale(IndexDoubleProcessor data, int dstart, int dend, double[] re, int rstart, int rend, int offset, UnivariateStatistic stat) {
		
		DoubleArrayList accumm = new DoubleArrayList();
		
		double factor = (dend-dstart)/(double)(rend-rstart);
		if (factor>=1) {
			int startInData = dstart;
			for (int p=rstart; p<Math.min(re.length-offset,rend); p++) {
				int endInData = dstart+(int) ((p-rstart+1)*factor);
				data.iterate((pos,val)->accumm.add(val),startInData,endInData);
//				re[p]/=endInData-startInData;
				accumm.add(0,endInData-startInData-accumm.size()); // this is important when stat=Mean (as the mean should be computed over all elements, not only over the non-zeros!)
				re[p+offset]+=stat.evaluate(accumm.getRaw(), 0, accumm.size());
				accumm.clear();
				startInData = endInData;
			}
		}
		else {
			data.iterate((pos,val)->{
				int s = (int) ((pos-dstart)/factor);
				int e = (int) ((pos+1-dstart)/factor);
				for (int p=s; p<e && rstart+p+offset<re.length; p++)
					re[rstart+p+offset]+=val;
			},dstart,dend);
		}
	}
	
	
	

	public static MetageneDefinition define(int width) {
		return new MetageneDefinition(width);
	}
	
	private static class MetageneRange {
		private String name;
		private int width;
		private MetageneRegionProvider reg;
		private MetageneDataProvider provider;
		private double where = -1;
		private int maxLength = -1;
		private double fixedFactor = -1;
		private UnivariateStatistic shrinkStat;
		
		public MetageneRange(int width,MetageneDefinition.MetageneRangeDefinition def) {
			this.width = width;
			this.name = def.name;
			this.reg = def.reg;
			this.provider = def.provider;
			this.where = def.where;
			this.maxLength = def.maxLength;
			this.fixedFactor = def.fixedFactor;
			this.shrinkStat = def.shrinkStat;
		}

	}

	
	public static void main(String[] args) throws IOException {
		
		Metagene meta = Metagene.define(600)
				.region("r1",MetageneDataProvider.getTestV())
					.setTruncateOrPad(0, 200)
					.setRegionProvider(MetageneRegionProvider.upstream(250))
					.add()
				.region("r2",MetageneDataProvider.getTestV())
					.setWidth(1)
					.setRegionProvider(MetageneRegionProvider.upstream(200))
					.add()
				.create();		
		
		ImmutableReferenceGenomicRegion<Transcript> t = null;
		
		meta.add("Test",ImmutableReferenceGenomicRegion.parse("1+:100-200"));
		
		meta.writeRds("test.rds");
		meta.writeTsv("test.tsv.gz");
	
	}


}
