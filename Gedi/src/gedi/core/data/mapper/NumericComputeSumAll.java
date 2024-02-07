package gedi.core.data.mapper;

import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.functions.NumericArrayFunction;

@GenomicRegionDataMapping(fromType=PixelBlockToValuesMap.class,toType=PixelBlockToValuesMap.class)
public class NumericComputeSumAll extends NumericCompute {

	public NumericComputeSumAll() {
		super(d->NumericArray.wrap(d.evaluate(NumericArrayFunction.Sum)));
	}

	
}
