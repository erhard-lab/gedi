package gedi.core.data.mapper;

import java.util.Collection;
import java.util.function.DoubleBinaryOperator;
import java.util.function.Function;
import java.util.function.UnaryOperator;

import javax.script.ScriptException;

import gedi.core.data.numeric.DenseGenomicNumericProvider;
import gedi.core.data.numeric.GenomicNumericProvider;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTuple;
import gedi.util.nashorn.JSFunction;

@GenomicRegionDataMapping(fromType=PixelBlockToValuesMap.class,toType=PixelBlockToValuesMap.class)
public class NumericCompute implements GenomicRegionDataMapper<PixelBlockToValuesMap, PixelBlockToValuesMap>{

	
	private Function<NumericArray,NumericArray> computer = t->t;
	
	
	public NumericCompute(String js) throws ScriptException {
		this.computer = new JSFunction(false, "function(data) "+js);
	}
	
	public NumericCompute(UnaryOperator<NumericArray> computer) {
		this.computer = computer;
	}

	@Override
	public PixelBlockToValuesMap map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			PixelBlockToValuesMap data) {
		
		
		NumericArray[] values = new NumericArray[data.size()];
		for (int i=0; i<values.length; i++)
			values[i] = computer.apply(data.getValues(i));
		
		
		return new PixelBlockToValuesMap(data, values);
		
	}


	
}
