package executables;

import gedi.lfc.localtest.javapipeline.LTestRunner;
import gedi.lfc.localtest.javapipeline.LTestParameterSet;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class LocalTest {

	
	public static void main(String[] args) {
		
	
		LTestParameterSet params = new LTestParameterSet();
		GediProgram pipeline = GediProgram.create("LocalTest",
				new LTestRunner(params)
				);
		GediProgram.run(pipeline, params.paramFile, new CommandLineHandler("LocalTest","LocalTest identifies regions in a gene that are differently regulated than all other regions.",args));
		
	}
}
