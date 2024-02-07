package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;

import gedi.grand3.experiment.ExperimentalDesign;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3WriteAllPseudobulkFile extends GediProgram {

	

	
	
	public Grand3WriteAllPseudobulkFile(Grand3ParameterSet params) {
		addInput(params.experimentalDesignFile);
		addOutput(params.allPseudobulkFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		File expDesign = getParameter(0);
		ExperimentalDesign design = ExperimentalDesign.fromTable(expDesign);

		LineWriter out = getOutputWriter(0);
		out.writeLine("Pseudobulk\tCell");
		for (int i=0; i<design.getCount(); i++) {
			out.writef("%s\t%s\n", design.getSampleName(i),design.getFullName(i));
		}
		out.close();
		return null;
	}


	
}
