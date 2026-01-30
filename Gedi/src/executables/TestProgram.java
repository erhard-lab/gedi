package executables;


import java.io.File;
import java.io.IOException;

import gedi.app.Gedi;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.javapipeline.Grand3BurstMcmcOutput;
import gedi.grand3.javapipeline.Grand3CollectSufficientStatistics;
import gedi.grand3.javapipeline.Grand3EstimateModel;
import gedi.grand3.javapipeline.Grand3OutputFlatfiles;
import gedi.grand3.javapipeline.Grand3ParameterSet;
import gedi.grand3.javapipeline.Grand3ProcessTargets;
import gedi.grand3.javapipeline.Grand3Resimulate;
import gedi.grand3.javapipeline.Grand3SetupTargetsAcceptor;
import gedi.grand3.javapipeline.Grand3SetupTargetsExonIntron;
import gedi.grand3.javapipeline.Grand3SetupTargetsGenes;
import gedi.grand3.javapipeline.Grand3SetupTargetsTest;
import gedi.grand3.javapipeline.Grand3SnpsAndClip;
import gedi.grand3.javapipeline.Grand3WriteAllPseudobulkFile;
import gedi.grand3.javapipeline.Grand3WriteExperimentalDesign;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.LogUtils.LogMode;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class TestProgram {

	public static void main(String[] args) {
		Gedi.startup(true,LogMode.Normal,"Test");
		
		TestParameterSet params = new TestParameterSet();
		
		GediProgram pipeline = GediProgram.create("Test",
				new TestFileReaderWriterProgram("init",params.init,params.prefix),
				new TestFileReaderWriterProgram("bA",params.brancha,params.prefix,params.init),
				new TestFileReaderWriterProgram("bB1",params.branchb1,params.prefix,params.init),
				new TestFileReaderWriterProgram("bB2",params.branchb2,params.prefix,params.branchb1),
				new TestFileReaderWriterProgram("result",params.result,params.prefix,params.brancha,params.branchb2)
				);
		
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("Test","Just a test for the program step system.",args));
	}

	
	public static class TestParameterSet extends GediParameterSet {

		public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call Grand3", false, new FileParameterType());
		public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());
		
		public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
		public GediParameter<File> init = new GediParameter<File>(this,"${prefix}.init", "", false, new FileParameterType());
		public GediParameter<File> brancha = new GediParameter<File>(this,"${prefix}.branchA", "", false, new FileParameterType());
		public GediParameter<File> branchb1 = new GediParameter<File>(this,"${prefix}.branchB1", "", false, new FileParameterType());
		public GediParameter<File> branchb2 = new GediParameter<File>(this,"${prefix}.branchB2", "", false, new FileParameterType());
		public GediParameter<File> result = new GediParameter<File>(this,"${prefix}.result", "", false, new FileParameterType());

	}
	
	
	public static class TestFileReaderWriterProgram extends GediProgram {

		private String add;

		public TestFileReaderWriterProgram(String add, GediParameter<File> output, GediParameter<String> prefix, GediParameter<File>...inputs) {
			super(add);
			this.add=add;
			this.addInput(prefix);
			for (GediParameter<File> i : inputs)
				addInput(i);
			addOutput(output);
		}
		
		
		
		public String execute(GediProgramContext context) throws IOException, InterruptedException {
			StringBuilder sb = new StringBuilder();
			for (int i=1; i<getInputSpec().size(); i++)
				sb.append(FileUtils.readAllText(getParameter(i)));
			sb.append(add);
			
			FileUtils.writeAllText(sb.toString(), getOutputFile(0));
			
			return null;
		}


		
	}
}
