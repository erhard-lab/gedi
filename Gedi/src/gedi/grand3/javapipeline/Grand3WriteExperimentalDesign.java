package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.util.FileUtils;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3WriteExperimentalDesign extends GediProgram {

	
	
	public Grand3WriteExperimentalDesign(Grand3ParameterSet params) {
		addInput(params.reads);
		addInput(params.auto);
		addOutput(params.experimentalDesignFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(pind++);
		boolean auto = getBooleanParameter(pind++);
		
		if (reads.getRandomRecord() instanceof BarcodedAlignedReadsData) {
			String bcFile = FileUtils.getExtensionSibling(reads.getPath(), "barcodes.tsv");
			ExperimentalDesign design = ExperimentalDesign.infer(context.getLog(),reads.getMetaDataConditions(), new File(bcFile));
			design.writeTable(getOutputFile(0));
		}
		else {
			ExperimentalDesign design = ExperimentalDesign.infer(context.getLog(),reads.getMetaDataConditions(),null);
			design.writeTable(getOutputFile(0));
		}
		
		// this is cool: Always exit here (except for -auto mode); this class is not executed if the file already exists!
		if (!auto) {
			System.out.println("Stopped analysis because of missing experimental design. Either restart (probably after editing the file "+getOutputFile(0)+" or use the -auto parameter to continue analysis!");
			System.exit(2);
		}
		
		return null;
	}


	
}
