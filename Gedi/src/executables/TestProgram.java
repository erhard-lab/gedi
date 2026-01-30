package executables;


import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

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
				new TestFileReaderWriterProgram("init",params.prefix,Arrays.asList(),Arrays.asList(params.init)),
				new TestFileReaderWriterProgram("bA",params.prefix,Arrays.asList(params.init),Arrays.asList(params.brancha,params.res2)),
				new TestFileReaderWriterProgram("bB1",params.prefix,Arrays.asList(params.init),Arrays.asList(params.branchb1)),
				new TestFileReaderWriterProgram("bB2",params.prefix,Arrays.asList(params.branchb1),Arrays.asList(params.branchb2)),
				new TestFileReaderWriterProgram("result",params.prefix,Arrays.asList(params.brancha,params.branchb2),Arrays.asList(params.result))
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
		public GediParameter<File> res2 = new GediParameter<File>(this,"${prefix}.res2", "", false, new FileParameterType());
		public GediParameter<File> result = new GediParameter<File>(this,"${prefix}.result", "", false, new FileParameterType());

	}
	
	
	public static class TestFileReaderWriterProgram extends GediProgram {

		private String add;

		public TestFileReaderWriterProgram(String add, GediParameter<String> prefix, List<GediParameter<File>> inputs, List<GediParameter<File>> outputs) {
			super(add);
			this.add=add;
			this.addInput(prefix);
			for (GediParameter<File> i : inputs)
				addInput(i);
			
			for (GediParameter<File> o : outputs)
				addOutput(o);
		}
		
		
		
		public String execute(GediProgramContext context) throws IOException, InterruptedException {
			StringBuilder sb = new StringBuilder();
			for (int i=1; i<getInputSpec().size(); i++)
				sb.append(FileUtils.readAllText(getParameter(i)));
			sb.append(add);
			
			for (int i=0; i<getOutputSpec().size(); i++)
				FileUtils.writeAllText(sb.toString(), getOutputFile(i));
			
			return null;
		}


		
	}
}
