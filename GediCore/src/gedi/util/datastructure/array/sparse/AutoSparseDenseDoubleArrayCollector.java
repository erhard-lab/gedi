/**
 * 
 *    Copyright 2017 Florian Erhard
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 * 
 */
package gedi.util.datastructure.array.sparse;

import java.io.IOException;
import java.util.Arrays;

import gedi.app.Config;
import gedi.util.ArrayUtils;
import gedi.util.functions.IntDoubleToDoubleFunction;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class AutoSparseDenseDoubleArrayCollector implements BinarySerializable {
	
	private static final double MARGIN = 2;
	
	private int length;
	
	private int[] indices;
	private double[] values;
	private int count;
	private boolean aggregated = true;
	
	private double[] dense;

	public AutoSparseDenseDoubleArrayCollector() {}
	
	public AutoSparseDenseDoubleArrayCollector(int expectedCapacity, int length) {
		int initialCap = (int) (expectedCapacity*MARGIN);
		this.length = length;
		
		if (length*Double.BYTES<initialCap*(Integer.BYTES+Double.BYTES)) {
			// start with dense!
			dense = new double[length];
		}
		else {
			indices = new int[initialCap];
			values = new double[initialCap];
		}
	}
	
	/**
	 * Iterates over all non-zero entries (ordered) and applies op.
	 * @param op
	 */
	public void process(IntDoubleToDoubleFunction op) {
		if (dense!=null) {
			for (int i=0; i<dense.length; i++)
				if (dense[i]!=0)
					dense[i] = op.applyAsDouble(i, dense[i]);
		}
		else {
			aggregate();
			for (int i=0; i<count; i++) 
				values[i] = op.applyAsDouble(indices[i], values[i]);
		}
	}
	
	
	public double get(int index) {
		if (dense!=null)
			return dense[index];
		aggregate();
		int ind = Arrays.binarySearch(indices, 0,count,index);
		if (ind<0) return 0;
		return values[ind];
	}
	
	public int length() {
		return length;
	}
	
	public void add(AutoSparseDenseDoubleArrayCollector o) {
		if (this.dense!=null && o.dense!=null) {
			for (int i=0; i<dense.length; i++)
				this.dense[i]+=o.dense[i];
		}
		else if (this.dense!=null) {
			for (int i=0; i<o.indices.length; i++) 
				this.dense[o.indices[i]]+=o.values[i];
		}
		else if (o.dense!=null) {
			this.dense = o.dense.clone();
			for (int i=0; i<indices.length; i++) 
				this.dense[indices[i]]+=values[i];
			this.indices = null;
			this.values = null;
		}
		else {
			for (int i=0; i<o.indices.length; i++) 
				add(o.indices[i],o.values[i]);
		}
	}
	
	
	public void add(int index, double value) {
		if (index<0 || index>=length) throw new ArrayIndexOutOfBoundsException(index);
		if (value==0) return;
		
		if (dense!=null) {
			dense[index]+=value;
		}
		else {
			if (count==indices.length) {
				
				aggregate();
				
				if (count==indices.length) {
					int newCapacity = indices.length + (indices.length >> 1) + 1;
					if (length*Double.BYTES<newCapacity*(Integer.BYTES+Double.BYTES)) {
						// switch to dense!
						enforceDense();
						dense[index]+=value;
					}
					else {
						int[] nindices = new int[newCapacity];
						double[] nvalues = new double[newCapacity];
						System.arraycopy(indices, 0, nindices, 0, count);
						System.arraycopy(values, 0, nvalues, 0, count);
						this.indices = nindices;
						this.values = nvalues;
						if (count>0 && indices[count-1]>=index)
							aggregated = false;
						indices[count]=index;
						values[count++]=value;
					}
					
				}
				else {
					if (count>0 && indices[count-1]>=index)
						aggregated = false;
					indices[count]=index;
					values[count++]=value;
				}
			}
			else {
				if (count>0 && indices[count-1]>=index)
					aggregated = false;
				indices[count]=index;
				values[count++]=value;
			}
		}
	}

	public double[] enforceDense() {
		dense = getDense();
		this.indices = null;
		this.values = null;
		return dense;
	}

	private double[] getDense() {
		if (dense!=null) return dense;
		
		double[] re = new double[length];
		for (int i=0; i<indices.length; i++) 
			re[indices[i]]+=values[i];
		
		return re;
	}


	private void aggregate() {
		if (dense==null && !aggregated) {
			
			ArrayUtils.parallelSort(indices, values, 0, count);
			int index = 0;
			for (int i=index+1; i<count; i++) {
				if (indices[i]!=indices[index]) {
					index++;
					values[index]=values[i];
					indices[index]=indices[i];
				} else 
					values[index]+=values[i];
			}
			count = index+1;
			
			aggregated = true;
		}
	}

	
	@Override
	public String toString() {
		aggregate();
		StringBuilder sb = new StringBuilder();
		sb.append("[");
		if (dense!=null) {
			for (int i=0; i<dense.length; i++) {
				if (i>0) sb.append(",");
				sb.append(String.format(Config.getInstance().getRealFormat(),dense[i]));
			}
		}
		else {
			for (int i=0; i<count; i++) {
				if (i>0) sb.append(",");
				sb.append(indices[i]).append(":").append(String.format(Config.getInstance().getRealFormat(),values[i]));
			}
		}
		sb.append("]");
		return sb.toString();
	}
	
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		aggregate();
		out.putInt(dense==null?length:-length);
		if (dense==null) {
			out.putCInt(count);
			for (int i=0; i<count; i++)
				out.putCInt(indices[i]);
			for (int i=0; i<count; i++)
				out.putDouble(values[i]);
		}
		else {
			for (int i=0; i<dense.length; i++)
				out.putDouble(dense[i]);
		}
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		dense = null;
		this.length = in.getInt();
		if (this.length<0) {
			this.length = -this.length;
			// dense mode
			this.dense = new double[length];
			for (int i=0; i<length; i++)
				this.dense[i] = in.getDouble();
		}
		else {
			count = in.getCInt();
			this.indices = new int[count];
			this.values = new double[count];
			for (int i=0; i<count; i++)
				indices[i]=in.getCInt();
			for (int i=0; i<count; i++)
				values[i]=in.getDouble();
		}
	}

	public static AutoSparseDenseDoubleArrayCollector wrap(double... a) {
		AutoSparseDenseDoubleArrayCollector re = new AutoSparseDenseDoubleArrayCollector(Math.max(50, a.length>>3),a.length);
		for (int i=0; i<a.length; i++) 
			re.add(i, a[i]);
		return re;
	}
	
}
