package gedi.slam;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Logger;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.MemoryDoubleArray;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.sparse.AutoSparseDenseIntArrayCollector;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class GeneData implements BinarySerializable {

	private String gene;
	private ReadData[] reads;
	private ReadData[] doublereads;
	private ReadData[] both;
	private NumericArray readCount;

	public GeneData(String gene, ReadData[] reads) {
		this.gene = gene;
		this.reads = reads;
		this.doublereads = new ReadData[0];
		this.readCount = NumericArray.createMemory(reads[0].getCount().length(), NumericArrayType.Double);
	}
	
	public GeneData(String gene, ReadData[] reads, NumericArray readCount) {
		this.gene = gene;
		this.reads = reads;
		this.doublereads = new ReadData[0];
		this.readCount = readCount;
	}
	
	public GeneData(String gene, ReadData[] reads, ReadData[] doublereads, ReadData[] both,NumericArray readCount) {
		this.gene = gene;
		this.reads = reads;
		this.doublereads = doublereads;
		this.both = both;
		this.readCount = readCount;
	}

	public GeneData() {
	}

	@Override
	public String toString() {
		return gene+"\t"+readCount.toArrayString(",", false);
	}

	public String getGene() {
		return gene;
	}

	public ReadData[] getReads() {
		return reads;
	}

	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putString(gene);
		readCount.serialize(out);
		out.putCInt(reads.length);
		for (int i=0; i<reads.length; i++)
			reads[i].serialize(out);
		out.putCInt(doublereads.length);
		for (int i=0; i<doublereads.length; i++)
			doublereads[i].serialize(out);
		out.putCInt(both.length);
		for (int i=0; i<both.length; i++)
			both[i].serialize(out);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		gene = in.getString();
		readCount = new MemoryDoubleArray();
		readCount.deserialize(in);
		reads = new ReadData[in.getCInt()];
		for (int i=0; i<reads.length; i++){
			reads[i] = new ReadData();
			reads[i].deserialize(in);
		}
		doublereads = new ReadData[in.getCInt()];
		for (int i=0; i<doublereads.length; i++){
			doublereads[i] = new ReadData();
			doublereads[i].deserialize(in);
		}
		both = new ReadData[in.getCInt()];
		for (int i=0; i<both.length; i++){
			both[i] = new ReadData();
			both[i].deserialize(in);
		}
	}
	
	
	private static int[] getConversionVector(int cond, ReadData[] reads) {
		int[] re = new int[reads.length];
		int index = 0;
		for (ReadData rd : reads)
			if (rd.getCount().get(cond)>0)
				re[index++] = rd.getConversions();
		return ArrayUtils.redimPreserve(re, index);
	}
	
	private static int[] getTotalVector(int cond, ReadData[] reads) {
		int[] re = new int[reads.length];
		int index = 0;
		for (ReadData rd : reads)
			if (rd.getCount().get(cond)>0)
				re[index++] = rd.getTotal();
		return ArrayUtils.redimPreserve(re, index);
	}
	
	private static double[] getWeightVector(int cond, ReadData[] reads) {
		double[] re = new double[reads.length];
		int index = 0;
		for (ReadData rd : reads)
			if (rd.getCount().get(cond)>0)
				re[index++] = rd.getCount().get(cond);
		return ArrayUtils.redimPreserve(re, index);
	}
	
	
	public int[] getConversionVector(int cond) {
		return getConversionVector(cond, reads);
	}
	
	public int[] getTotalVector(int cond) {
		return getTotalVector(cond,reads);
	}
	
	public double[] getWeightVector(int cond) {
		return getWeightVector(cond, reads);
	}
	
	public int[] getConversionDoubleVector(int cond) {
		return getConversionVector(cond, doublereads);
	}
	
	public int[] getTotalDoubleVector(int cond) {
		return getTotalVector(cond,doublereads);
	}
	
	public double[] getWeightDoubleVector(int cond) {
		return getWeightVector(cond, doublereads);
	}
	
	public int[] getConversionBothVector(int cond) {
		return getConversionVector(cond, both);
	}
	
	public int[] getTotalBothVector(int cond) {
		return getTotalVector(cond,both);
	}
	
	public double[] getWeightBothVector(int cond) {
		return getWeightVector(cond, both);
	}

	
	private int dim() {
		if (reads.length>0) return reads[0].getCount().length();
		if (doublereads.length>0) return doublereads[0].getCount().length();
		return readCount.length();
	}
	
	public NumericArray getReadCount() {
		return readCount;
	}
	
	
	public NumericArray getTotalConversions() {
		NumericArray re = NumericArray.createMemory(dim(),NumericArrayType.Double);
		for (int i=0; i<re.length(); i++) {
			int cond = i;
			double count =  EI.wrap(reads).mapToDouble(r->r.getConversions()*r.getCount().get(cond)).sum();
			re.setDouble(i, count);
		}
		return re;
	}
	
	public NumericArray getTotalCoverage() {
		NumericArray re = NumericArray.createMemory(dim(),NumericArrayType.Double);
		for (int i=0; i<re.length(); i++) {
			int cond = i;
			double total =  EI.wrap(reads).mapToDouble(r->r.getTotal()*r.getCount().get(cond)).sum();
			re.setDouble(i, total);
		}
		return re;
	}
	
	public NumericArray getTotalDoubleHits() {
		NumericArray re = NumericArray.createMemory(dim(),NumericArrayType.Double);
		for (int i=0; i<re.length(); i++) {
			int cond = i;
			double count =  EI.wrap(doublereads).mapToDouble(r->r.getConversions()*r.getCount().get(cond)).sum();
			re.setDouble(i, count);
		}
		return re;
	}
	
	public NumericArray getTotalDoubleHitCoverage() {
		NumericArray re = NumericArray.createMemory(dim(),NumericArrayType.Double);
		for (int i=0; i<re.length(); i++) {
			int cond = i;
			double total =  EI.wrap(doublereads).mapToDouble(r->r.getTotal()*r.getCount().get(cond)).sum();
			re.setDouble(i, total);
		}
		return re;
	}

	public NumericArray getReadsWithConversions(int minconv) {
		return EI.wrap(reads).filter(r->r.getConversions()>=minconv).reduce(NumericArray.createMemory(dim(),NumericArrayType.Double), (n,r)->{
			r.add(n.getCount());
			return r;
		});
	}

	public GeneProportion infer(OptimNumericalIntegrationProportion[] vb, Logger log) {
		int cond = vb.length;
		
		// number of single and double hits
		int[] singleHits = new int[cond];
		int[] doubleHits = new int[cond];
		
		
		for (ReadData rd : reads) 
			rd.getCount().process((index,value)->{if (value>0) singleHits[index]++; return value;});
		for (ReadData rd : doublereads) 
			rd.getCount().process((index,value)->{if (value>0) doubleHits[index]++; return value;});
		
		int[][] d = new int[cond][];
		int[][] n = new int[cond][];
		double[][] m = new double[cond][];
		for (int i=0; i<cond; i++) {
			d[i] = new int[singleHits[i]];
			n[i] = new int[singleHits[i]];
			m[i] = new double[singleHits[i]];
		}
		
		int[][] d2 = new int[cond][];
		int[][] n2 = new int[cond][];
		double[][] m2 = new double[cond][];
		for (int i=0; i<cond; i++) {
			d2[i] = new int[doubleHits[i]];
			n2[i] = new int[doubleHits[i]];
			m2[i] = new double[doubleHits[i]];
		}
		
		Arrays.fill(singleHits, 0);
		Arrays.fill(doubleHits, 0);
		
		for (ReadData rd : reads) {
			rd.getCount().process((index,value)->{
				d[index][singleHits[index]]=rd.getConversions();
				n[index][singleHits[index]]=rd.getTotal();
				m[index][singleHits[index]]=value;
				singleHits[index]++;
				return value;
			});
		}
			
		for (ReadData rd : doublereads) {
			rd.getCount().process((index,value)->{
				d2[index][doubleHits[index]]=rd.getConversions();
				n2[index][doubleHits[index]]=rd.getTotal();
				m2[index][doubleHits[index]]=value;
				doubleHits[index]++;
				return value;
			});
		}
			
		SlamEstimationResult[] re = new SlamEstimationResult[vb.length];
		for (int i=0; i<cond; i++)
			if (vb[i]!=null)
				re[i] = vb[i].infer(d[i], n[i], m[i], d2[i], n2[i], m2[i], log);
		
		return new GeneProportion(this, re);
	}
	
	
}
