package gedi.riboseq.codonprocessor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Consumer;
import java.util.function.Function;

import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.mutable.MutableTuple;

public class ExpandingCodonProcessorCounter extends CodonProcessorCounter {

	
	
	
	private Function<HashMap<String, Object>, ExtendedIterator<HashMap<String,Object>>> expander;
	
	public ExpandingCodonProcessorCounter(String prefix, Function<HashMap<String,Object>,ExtendedIterator<HashMap<String,Object>>> expander, String... varNames)  {
		super(prefix,varNames);
		this.expander = expander;
	}

	@Override
	public void count(HashMap<String, Object> vars, NumericArray count) {
		
		for (HashMap<String,Object> sub : expander.apply(vars).loop()) {
			if (testConditions(sub)) {
				setup(sub);
				count(count);
			}
		}
		
	}


	
		
	
}
