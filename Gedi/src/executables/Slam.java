package executables;


import gedi.app.Gedi;
import gedi.slam.javapipeline.SlamCollectStatistics;
import gedi.slam.javapipeline.SlamDetectSnps;
import gedi.slam.javapipeline.SlamInfer;
import gedi.slam.javapipeline.SlamParameterSet;
import gedi.slam.javapipeline.SlamRates;
import gedi.slam.javapipeline.SlamTest;
import gedi.util.LogUtils.LogMode;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class Slam {



	private static String getChangelog() {
		return "2.0.6:\n"
				+ "introduced the -modelall parameter\n"
				+ "added plotting of ntr for overall read types\n"
				+ "added ${prefix}.strandness file for easier restarts\n"
				+ "introduced -highmem parameter\n\n"
				+ "2.0.7: geneData/extData/snpData files are removed in the end, tsv is gzipped\n\n";
	}

	
	public static void main(String[] args) {
		Gedi.startup(true,LogMode.Normal,"GRAND-SLAM");
		
		SlamParameterSet params = new SlamParameterSet();
		GediProgram pipeline = GediProgram.create("GrandSlam",
				new SlamDetectSnps(params),
				new SlamCollectStatistics(params),
				new SlamRates(params),
				new SlamInfer(params),
				new SlamTest(params)
				);
		pipeline.setChangelog(getChangelog());
		
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("GrandSlam","Grand-Slam is an analysis method for SLAM/scSLAM-seq data.",args));
		
	}
	
}
