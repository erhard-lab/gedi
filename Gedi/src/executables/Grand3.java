package executables;


import gedi.app.Gedi;
import gedi.grand3.javapipeline.Grand3BurstMcmcOutput;
import gedi.grand3.javapipeline.Grand3CollectSufficientStatistics;
import gedi.grand3.javapipeline.Grand3EstimateModel;
import gedi.grand3.javapipeline.Grand3OutputFlatfiles;
import gedi.grand3.javapipeline.Grand3ParameterSet;
import gedi.grand3.javapipeline.Grand3ProcessTargets;
import gedi.grand3.javapipeline.Grand3Resimulate;
import gedi.grand3.javapipeline.Grand3SetupTargetsAcceptor;
import gedi.grand3.javapipeline.Grand3SetupTargetsGenes;
import gedi.grand3.javapipeline.Grand3SetupTargetsTest;
import gedi.grand3.javapipeline.Grand3SnpsAndClip;
import gedi.grand3.javapipeline.Grand3WriteAllPseudobulkFile;
import gedi.grand3.javapipeline.Grand3WriteExperimentalDesign;
import gedi.util.ArrayUtils;
import gedi.util.LogUtils.LogMode;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class Grand3 {

	public static void main(String[] args) {
		Gedi.startup(true,LogMode.Normal,"GRAND3");
		
		Grand3ParameterSet params = new Grand3ParameterSet();
		boolean hasMappedTarget = ArrayUtils.find(args, "-"+params.pseudobulkFile.getName())>=0;
		boolean hasTargetMixmat = ArrayUtils.find(args, "-"+params.targetMixmat.getName())>=0;
		
		GediProgram pipeline = GediProgram.create("Grand3",
				new Grand3WriteExperimentalDesign(params),
				new Grand3SetupTargetsAcceptor(params),
				new Grand3SetupTargetsGenes(params),
				new Grand3SetupTargetsTest(params),
				new Grand3SnpsAndClip(params),
				new Grand3CollectSufficientStatistics(params),
				new Grand3EstimateModel(params),
				new Grand3ProcessTargets(params, hasMappedTarget, hasTargetMixmat),
				hasTargetMixmat? null: new Grand3OutputFlatfiles(params, hasMappedTarget),
				new Grand3WriteAllPseudobulkFile(params),
				hasTargetMixmat? null: new Grand3Resimulate<>(params, hasMappedTarget),
				new Grand3BurstMcmcOutput<>(params)
				);
		
		pipeline.setChangelog(getChangelog());
		
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("Grand3","Grand3 is an analysis method for (sc-) SLAM/TimeLapse/TUC-seq data.",args));
	}

	private static String getChangelog() {
		return "3.0.3:\n"
				+ " introduced the -restrict-subreads parameter\n\n"+
				"3.0.4:\n"
				+ " implemented 5' sequencing (10x): antisense with UMIs\n per gene output for the MixMatrices\n\n"+
				"3.0.5:\n"
				+ " implemented output for splice acceptor sites\n\n"
				;
	}
		
	
}
