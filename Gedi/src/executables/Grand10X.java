package executables;


import java.io.IOException;

import gedi.grand10x.javapipeline.Grand10XClusterProgram;
import gedi.grand10x.javapipeline.Grand10XDemultiplexProgram;
import gedi.grand10x.javapipeline.Grand10xClusterPositionInferenceProgram;
import gedi.grand10x.javapipeline.Grand10xParameterSet;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;


public class Grand10X {
	
	




	public static void main(String[] args) throws IOException {
		
		
		Grand10xParameterSet params = new Grand10xParameterSet();
		GediProgram pipeline = GediProgram.create("Grand10X",
				new Grand10xClusterPositionInferenceProgram(params),
				new Grand10XClusterProgram(params),
				new Grand10XDemultiplexProgram(params)
				);
		GediProgram.run(pipeline, null, null, new CommandLineHandler("Grand10X","Grand10X reads 10x cit files and creates a demultiplexed and de-barcoded cit file.",args));

	}
	
}
