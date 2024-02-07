package executables;


import gedi.app.Gedi;
import gedi.grand3.javapipeline.Grand3EstimateModel;
import gedi.grand3.javapipeline.Grand3ParameterSet;
import gedi.util.ArrayUtils;
import gedi.util.LogUtils.LogMode;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class Grand3Test {

	public static void main(String[] args) {
		Gedi.startup(true,LogMode.Normal,"GRAND3");
		
		Grand3ParameterSet params = new Grand3ParameterSet();
		boolean hasMappedTarget = ArrayUtils.find(args, "-"+params.pseudobulkFile.getName())>=0;
		
		GediProgram pipeline = GediProgram.create("Grand3",
				new Grand3EstimateModel(params)
				);
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("Grand3","Grand3 is an analysis method for (sc-) SLAM/TimeLapse/TUC-seq data.",args));
	}
		
	
}
