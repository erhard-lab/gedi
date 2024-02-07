package gedi.util.datastructure.array.computed;

import java.io.IOException;
import java.util.function.IntToDoubleFunction;
import java.util.function.IntUnaryOperator;

import gedi.util.datastructure.array.DoubleArray;
import gedi.util.datastructure.array.IntegerArray;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.io.randomaccess.BinaryReader;

public class ComputedDoubleArray extends DoubleArray {

	private IntToDoubleFunction function;
	private int length;
	
	public ComputedDoubleArray(IntToDoubleFunction function, int length) {
		this.function = function;
		this.length = length;
	}

	@Override
	public void setDouble(int index, double value) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isReadOnly() {
		return true;
	}

	@Override
	public int length() {
		return length;
	}

	@Override
	public double getDouble(int index) {
		return function.applyAsDouble(index);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public NumericArray clear() {
		throw new UnsupportedOperationException();
	}
	
}
