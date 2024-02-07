package gedi.grand3.knmatrix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.EI;
import gedi.util.functions.ParallelizedState;
import gedi.util.io.text.LineWriter;

public class KNMatrix implements ParallelizedState<KNMatrix> {

	private TreeMap<Integer,KVec> mat = new TreeMap<>();
	private int numCond;
	
	
	public KNMatrix(int numCond) {
		this.numCond = numCond;
	}
	
	public AutoSparseDenseDoubleArrayCollector get(int k, int n) {
		KVec a = mat.computeIfAbsent(n, x->new KVec(x)); 
		return a.get(k);
	}

	@Override
	public KNMatrix spawn(int index) {
		return new KNMatrix(numCond);
	}
	
	public int getNumConditions() {
		return numCond;
	}

	@Override
	public void integrate(KNMatrix other) {
		Iterator<Entry<Integer, KVec>> it = other.mat.entrySet().iterator();
		while (it.hasNext()) {
			Entry<Integer,KVec> en = it.next();
			mat.compute(en.getKey(), (n,present)->{
				if (present==null) return en.getValue();
				present.integrate(en.getValue());
				return present;
			});
		}
	}
	
	public void clear() {
		mat = new TreeMap<>();
	}

	
	public void write(String prefix, int condition, LineWriter wr) throws IOException {
		Iterator<Entry<Integer, KVec>> it = mat.entrySet().iterator();
		while (it.hasNext()) {
			Entry<Integer,KVec> en = it.next();
			Integer n = en.getKey();
			
			for (int k=0; k<en.getValue().densePart.length; k++) {
				double count = en.getValue().densePart[k].get(condition);
				if (count>0)
					wr.writef("%s%d\t%d\t%.0f\n", prefix,k,n,count);
			}
			
			Iterator<Entry<Integer, AutoSparseDenseDoubleArrayCollector>> kit = en.getValue().sparsePart.entrySet().iterator();
			while (kit.hasNext()) {
				Entry<Integer,AutoSparseDenseDoubleArrayCollector> ken = kit.next();
				Integer k = ken.getKey();
				double count = ken.getValue().get(condition);
				wr.writef("%s%d\t%d\t%.0f\n", prefix,k,n,count);
			}
			
		}
	}

	public double[] computeExpectedProbs() {
		double[][] all = EI.wrap(mat.values()).map(v->v.computeExpectedProbsFactors()).reduce((a,b)->{
			ArrayUtils.add(a[0], b[0]);
			ArrayUtils.add(a[1], b[1]);
			return a;
		});
		
		for (int c=0; c<all.length; c++)
			all[0][c]=all[0][c]/all[1][c];
		
		return all[0];
	}

	/**
	 * condition,n->count
	 * @return
	 */
	public int[][] computeNvectors() {
		return  EI.wrap(mat.values()).map(v->v.getTotals()).toArray(int[].class);
	}
	
	
	/**
	 * condition,list
	 * @return
	 */
	public KNMatrixElement[][] toElements() {
		ArrayList<KNMatrixElement>[] re = new ArrayList[numCond];
		for (int c=0; c<numCond; c++) 
			re[c] = new  ArrayList<>();
		
		Iterator<Entry<Integer, KVec>> it = mat.entrySet().iterator();
		while (it.hasNext()) {
			Entry<Integer,KVec> en = it.next();
			Integer n = en.getKey();
			
			for (int k=0; k<en.getValue().densePart.length; k++) {
				AutoSparseDenseDoubleArrayCollector cv = en.getValue().densePart[k];
				int uk = k;
				cv.process((cond,count)->{
					if (count>0)
						re[cond].add(new KNMatrixElement(uk,n,count));
					return count;
				});
			}
			
			Iterator<Entry<Integer, AutoSparseDenseDoubleArrayCollector>> kit = en.getValue().sparsePart.entrySet().iterator();
			while (kit.hasNext()) {
				Entry<Integer,AutoSparseDenseDoubleArrayCollector> ken = kit.next();
				Integer k = ken.getKey();
				AutoSparseDenseDoubleArrayCollector cv = ken.getValue();
				cv.process((cond,count)->{
					if (count>0)
						re[cond].add(new KNMatrixElement(k,n,count));
					return count;
				});
			}
			
		}
		
		KNMatrixElement[][] re2 = new KNMatrixElement[numCond][];
		for (int i=0; i<re.length; i++)
			re2[i] = re[i].toArray(new KNMatrixElement[0]);
		
		return re2;
	}
	
	
	public KNMatrixElement[] toElements(int condition) {
		ArrayList<KNMatrixElement> re = new ArrayList<>();
		
		Iterator<Entry<Integer, KVec>> it = mat.entrySet().iterator();
		while (it.hasNext()) {
			Entry<Integer,KVec> en = it.next();
			Integer n = en.getKey();
			
			for (int k=0; k<en.getValue().densePart.length; k++) {
				AutoSparseDenseDoubleArrayCollector cv = en.getValue().densePart[k];
				double count = cv.get(condition);
				if (count>0) re.add(new KNMatrixElement(k, n, count));
			}
			
			Iterator<Entry<Integer, AutoSparseDenseDoubleArrayCollector>> kit = en.getValue().sparsePart.entrySet().iterator();
			while (kit.hasNext()) {
				Entry<Integer,AutoSparseDenseDoubleArrayCollector> ken = kit.next();
				Integer k = ken.getKey();
				AutoSparseDenseDoubleArrayCollector cv = ken.getValue();
				double count = cv.get(condition);
				if (count>0) re.add(new KNMatrixElement(k, n, count));
			}
			
		}
		
		return re.toArray(new KNMatrixElement[re.size()]);
	}
	
	
	private class KVec {
		int n;
		AutoSparseDenseDoubleArrayCollector[] densePart;
		TreeMap<Integer, AutoSparseDenseDoubleArrayCollector> sparsePart;
		public KVec(int n) {
			this.n = n;
			densePart = new AutoSparseDenseDoubleArrayCollector[10];
			for (int i=0; i<densePart.length; i++)
				densePart[i] = new AutoSparseDenseDoubleArrayCollector(numCond<50?numCond:10, numCond);
			sparsePart = new TreeMap<Integer, AutoSparseDenseDoubleArrayCollector>();
		}
		public int[] getTotals() {
			double[] re = new double[densePart[0].length()];
			for (int k=0; k<densePart.length; k++) {
				for (int c=0; c<re.length; c++) {
					re[c]+=densePart[k].get(c);
				}
			}
			for (Integer k : sparsePart.keySet()) {
				AutoSparseDenseDoubleArrayCollector pp = sparsePart.get(k);
				for (int c=0; c<re.length; c++) {
					re[c]+=pp.get(c);
				}
			}
			return ArrayUtils.toInt(re);
		}
		
		public double[][] computeExpectedProbsFactors() {
			double[] num = new double[densePart[0].length()];
			double[] denom = new double[densePart[0].length()];
			for (int k=0; k<densePart.length; k++) {
				for (int c=0; c<num.length; c++) {
					num[c]+=k*densePart[k].get(c);
					denom[c]+=densePart[k].get(c);
				}
			}
			for (Integer k : sparsePart.keySet()) {
				AutoSparseDenseDoubleArrayCollector pp = sparsePart.get(k);
				for (int c=0; c<num.length; c++) {
					num[c]+=k*pp.get(c);
					denom[c]+=pp.get(c);
				}
			}
			ArrayUtils.mult(num, 1.0/n);
			return new double[][] {num,denom};
		}
		public void integrate(KVec other) {
			for (int i=0; i<densePart.length; i++)
				densePart[i].add(other.densePart[i]);
			
			Iterator<Entry<Integer, AutoSparseDenseDoubleArrayCollector>> it = other.sparsePart.entrySet().iterator();
			while (it.hasNext()) {
				Entry<Integer,AutoSparseDenseDoubleArrayCollector> en = it.next();
				sparsePart.compute(en.getKey(), (n,present)->{
					if (present==null) return en.getValue();
					present.add(en.getValue());
					return present;
				});
			}
		}
		public AutoSparseDenseDoubleArrayCollector get(int k) {
			return k<densePart.length?densePart[k]:sparsePart.computeIfAbsent(k, x->new AutoSparseDenseDoubleArrayCollector(numCond<50?numCond:10, numCond));
		}
	}


	public static class KNMatrixElement {
		public int k;
		public int n;
		public double count;
		public KNMatrixElement(int k, int n, double count) {
			this.k = k;
			this.n = n;
			this.count = count;
		}
		@Override
		public String toString() {
			return k+"|"+n+" "+String.format("%.0f", count);
		}
	}


	public static KNMatrixElement[] sum(KNMatrixElement[] a, KNMatrixElement[] b) {
		Comparator<KNMatrixElement> comp = (aa,bb)->{
			int re = Integer.compare(aa.n, bb.n);
			if (re==0) re = Integer.compare(aa.k, bb.k);
			return re;
		};
		// assert sorting
		if (!ArrayUtils.isStrictAscending(a, comp) || !ArrayUtils.isStrictAscending(b,comp)) throw new RuntimeException("Cannot be!");
		
		KNMatrixElement[] re = EI.wrap(a).merge(comp, EI.wrap(b)).multiplex(comp, KNMatrixElement.class)
			.map(arr->{
				if (arr.length==1) return arr[0];
				if (arr.length>2) throw new RuntimeException("Cannot be!");
				return new KNMatrixElement(arr[0].k, arr[0].n, arr[0].count+arr[1].count);
			}).toArray(KNMatrixElement.class);
		
		if (!ArrayUtils.isStrictAscending(re, comp) || Math.abs(EI.wrap(a).mapToDouble(aa->aa.count).sum()+EI.wrap(b).mapToDouble(aa->aa.count).sum()-EI.wrap(re).mapToDouble(aa->aa.count).sum())>0.01) throw new RuntimeException("Cannot be!");
		return re;
	}

	
	
}
