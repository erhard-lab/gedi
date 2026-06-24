package gedi.grand3.javapipeline;


import java.io.IOException;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediParameterSpec;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3LoadReadToMem<A extends AlignedReadsData> extends GediProgram {

	
	
	public Grand3LoadReadToMem(Grand3ParameterSet params) {
		addInput(params.reads);
		addInput(params.tomem);
		addOutput(params.reads);
		
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		return null;
	}


	@Override
	protected void initParameter(GediParameterSet params, GediParameterSpec inputSpec) {
		super.initParameter(params, inputSpec);
		
		if (getBooleanParameter(1)) {
			GenomicRegionStorage<A> reads = getParameter(0);
			setOutput(0, reads.toMemory());
		}
		
	}
}
