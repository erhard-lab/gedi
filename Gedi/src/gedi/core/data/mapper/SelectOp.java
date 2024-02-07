package gedi.core.data.mapper;

import java.util.function.UnaryOperator;

import gedi.util.ParseUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.collections.intcollections.IntArrayList;

class SelectOp implements UnaryOperator<NumericArray> {

	private String range;
	private int[] pos;

	public SelectOp(String range) {
		this.range = range;
	}

	@Override
	public NumericArray apply(NumericArray t) {
		if (pos==null) 
			pos = ParseUtils.parseRangePositions(range, t.length(), new IntArrayList()).toIntArray();
		
		NumericArray re = NumericArray.createMemory(pos.length, t.getType());
		for (int i=0; i<pos.length; i++)
			re.copy(t, pos[i], i);
		
		return re;
	}
	
}