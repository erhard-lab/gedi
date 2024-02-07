package gedi.riboseq.codonprocessor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.mutable.MutableTuple;

public class SimpleCodonProcessorCounter extends CodonProcessorCounter {

	
	
	public SimpleCodonProcessorCounter(String prefix, String...vars)  {
		super(prefix,vars);
	}
	

	
	@Override
	public void count(HashMap<String, Object> vars, NumericArray count) {
		if (testConditions(vars)) {
			setup(vars);
			count(count);
		}
	}

	
	
}
