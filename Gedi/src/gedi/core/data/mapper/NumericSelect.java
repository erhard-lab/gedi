package gedi.core.data.mapper;

import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;

@GenomicRegionDataMapping(fromType=PixelBlockToValuesMap.class,toType=PixelBlockToValuesMap.class)
public class NumericSelect extends NumericCompute {

	

	public NumericSelect(String range) {
		super(new SelectOp(range));
	}


	
}
