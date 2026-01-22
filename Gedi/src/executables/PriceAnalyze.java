package executables;

import gedi.app.Gedi;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.riboseq.javapipeline.analyze.PriceAnalyzeMajorIsoforms;
import gedi.riboseq.javapipeline.analyze.PriceGenerateCodonProfiles;
import gedi.riboseq.javapipeline.analyze.PriceIdentifyMajorIsoform;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class PriceAnalyze {

	
	public static void main(String[] args) {
		
		Gedi.startup(true);
		
		PriceParameterSet params = new PriceParameterSet();
		GediProgram pipeline = GediProgram.create("PriceAnalyze",
				new PriceIdentifyMajorIsoform(params),
				new PriceGenerateCodonProfiles(params),
				new PriceAnalyzeMajorIsoforms(params)
				);
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("PriceAnalyze","PriceAnalyze analyzes further aspects  in PRICE-inferred codons.",args));
		
	}
	
}
