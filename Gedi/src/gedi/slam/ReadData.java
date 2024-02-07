package gedi.slam;

import java.io.IOException;
import java.util.Arrays;

import gedi.util.FileUtils;
import gedi.util.datastructure.array.MemoryDoubleArray;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class ReadData implements BinarySerializable{
	private int total;
	private int conversions;
	private AutoSparseDenseDoubleArrayCollector count;
	
	public ReadData() {
	}
	public ReadData(int total, int conversions, AutoSparseDenseDoubleArrayCollector ds) {
		super();
		this.total = total;
		this.conversions = conversions;
		this.count = ds;
	}
	@Override
	public String toString() {
		return "ReadData [total=" + total + ", conversions=" + conversions + ", count=" + count + "]";
	}
	public int getTotal() {
		return total;
	}
	public int getConversions() {
		return conversions;
	}
	public AutoSparseDenseDoubleArrayCollector getCount() {
		return count;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + conversions;
		result = prime * result + count.hashCode();
		result = prime * result + total;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ReadData other = (ReadData) obj;
		if (conversions != other.conversions)
			return false;
		if (!count.equals(other.count))
			return false;
		if (total != other.total)
			return false;
		return true;
	}
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putCInt(total);
		out.putCInt(conversions);
		count.serialize(out);
	}
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		total = in.getCInt();
		conversions = in.getCInt();
		count = new AutoSparseDenseDoubleArrayCollector();
		count.deserialize(in);
	}
	public void add(ReadData readData) {
		count.add(readData.count);
	}
	
	
	
	
}
