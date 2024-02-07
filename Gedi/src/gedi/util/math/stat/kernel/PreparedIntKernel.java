package gedi.util.math.stat.kernel;

import gedi.util.ArrayUtils;
import gedi.util.functions.EI;

public class PreparedIntKernel implements Kernel {

	private Kernel kernel;
	private double[] weights;
	private double[] cumWeights;
	private int halfInt;
	
	public PreparedIntKernel(Kernel kernel, boolean normalize) {
		this.kernel = kernel;
		halfInt = (int)Math.floor(kernel.halfSize());
		weights = new double[halfInt*2+1];
		for (int i=0; i<weights.length; i++)
			weights[i] = kernel.applyAsDouble(i-halfInt);
		if (normalize) 
			ArrayUtils.normalize(weights);
		cumWeights = weights.clone();
		ArrayUtils.cumSumInPlace(cumWeights, 1);
	}
	
	
	@Override
	public String name() {
		return kernel.name();
	}
	
	@Override
	public String toString() {
		return name();
	}
	
	public int getMinAffectedIndex(int index) {
		return index-halfInt;
	}
	
	public int getMaxAffectedIndex(int index) {
		return index+halfInt;
	}
	
	@Override
	public double applyAsDouble(double operand) {
		if ((int)operand!=operand) throw new RuntimeException("Only ints allowed!");
		int p = (int)operand;
		if (p<-halfInt||p>halfInt) return 0;
		return weights[p+halfInt];
	}

	@Override
	public double halfSize() {
		return halfInt;
	}

	
	public double[] processParallelized(int nthreads, double[] a, int start, int end) {
		if (nthreads<=1 || end-start<=1+halfInt*2) {
			double[] re = a.clone();
			processInPlace(re, start, end);
			return re;
		}
		nthreads=Math.min(nthreads, (end-start)/halfInt);
		
		int[] bound = new int[nthreads+1];
		for (int i=1; i<=nthreads; i++)
			bound[i] = (int) ((end-start)*(double)i/nthreads);
		
		double[] re = new double[end-start];
		for (double[] part : EI.seq(0, nthreads).parallelized(nthreads, 1, ei->ei.map(i->{
						int s = bound[i];
						int e = bound[i+1];
						double[] re2 = new double[e-s+1+2*halfInt];
						if (s==0)
							System.arraycopy(a, s, re2, halfInt, e-s+halfInt);
						else if (e==re.length)
							System.arraycopy(a, s-halfInt, re2, 0, e-s+halfInt);
						else
							System.arraycopy(a, s-halfInt, re2, 0, e-s+2*halfInt);
						processInPlace(re2, halfInt, e-s+halfInt);
						re2[re2.length-1] = i;
						return re2;
					})).loop()) {
			int i = (int) part[part.length-1];
			int s = bound[i];
			int e = bound[i+1];
			System.arraycopy(part, halfInt, re, s, e-s);
		}
		return re;
	}
	
	public double applyToArraySlice(double[] a, int pos) {
		
		int start = Math.max(0, pos-halfInt);
		int stop = Math.min(a.length-1, pos+halfInt);
		
		double re = 0;
		for (int i=start; i<=stop; i++)
			re+=weights[i-pos+halfInt]*a[i];
		
		return re;
	}

	
	public double applyToArraySlice(int[] a, int pos) {
		
		int start = Math.max(0, pos-halfInt);
		int stop = Math.min(a.length-1, pos+halfInt);
		
		double re = 0;
		for (int i=start; i<=stop; i++)
			re+=weights[i-pos+halfInt]*a[i];
		
		return re;
	}
	public double applyToArraySlice(long[] a, int pos) {
		
		int start = Math.max(0, pos-halfInt);
		int stop = Math.min(a.length-1, pos+halfInt);
		
		double re = 0;
		for (int i=start; i<=stop; i++)
			re+=weights[i-pos+halfInt]*a[i];
		
		return re;
	}


	public void processInPlace(double[] a, int from, int to) {
		double[] buff = new double[halfInt+1];
		
		for (int i=from; i<to; i++) {
			buff[i%buff.length] = applyToArraySlice(a, i);
			if (i-halfInt>=from)
				a[i-halfInt] = buff[(i+1)%buff.length];
		}
		for (int i=to; i<to+Math.min(halfInt, to-from); i++)
			if (i-halfInt>=from)
				a[i-halfInt] = buff[(i+1)%buff.length];
	}
	
	
	public void processInPlace(int[] a, int from, int to, double fac) {
		int[] buff = new int[halfInt+1];
		
		for (int i=from; i<to; i++) {
			buff[i%buff.length] = (int) (fac*applyToArraySlice(a, i));
			if (i-halfInt>=from)
				a[i-halfInt] = buff[(i+1)%buff.length];
		}
		for (int i=to; i<to+Math.min(halfInt, to-from); i++)
			a[i-halfInt] = buff[(i+1)%buff.length];
	}

	public void processInPlace(long[] a, int from, int to, double fac) {
		long[] buff = new long[halfInt+1];
		
		for (int i=from; i<to; i++) {
			buff[i%buff.length] = (long) (fac*applyToArraySlice(a, i));
			if (i-halfInt>=from)
				a[i-halfInt] = buff[(i+1)%buff.length];
		}
		for (int i=to; i<to+Math.min(halfInt, to-from); i++)
			a[i-halfInt] = buff[(i+1)%buff.length];
	}

	
	
}
