package gedi.core.data.mapper;

import java.util.function.BiFunction;

import javax.script.ScriptException;

import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.gui.genovis.pixelMapping.PixelBlockToValuesMap;
import gedi.gui.genovis.pixelMapping.PixelLocationMapping;
import gedi.gui.genovis.pixelMapping.PixelLocationMappingBlock;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.mutable.MutableTuple;
import gedi.util.nashorn.JSBiFunction;

@GenomicRegionDataMapping(fromType=MutableTuple.class,toType=PixelBlockToValuesMap.class)
public class NumericMerge implements GenomicRegionDataMapper<MutableTuple, PixelBlockToValuesMap>{

	
	private BiFunction<NumericArray,PixelLocationMappingBlock,NumericArray> computer = null;
	
	
	public void setCompute(String js) throws ScriptException {
		this.computer = new JSBiFunction<>(false, "function(data,block) "+js);
	}
	

	@Override
	public PixelBlockToValuesMap map(ReferenceSequence reference,
			GenomicRegion region,PixelLocationMapping pixelMapping,
			MutableTuple data2) {
		PixelBlockToValuesMap data = new PixelBlockToValuesMap(data2);
		if (computer==null) 
			return data;
		
		NumericArray[] values = new NumericArray[data.size()];
		for (int i=0; i<values.length; i++)
			values[i] = computer.apply(data.getValues(i),data.getBlock(i));
		
		
		return new PixelBlockToValuesMap(data, values);
		
		
	}


	
}
