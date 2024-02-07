package executables;

import gedi.app.Gedi;
import gedi.riboseq.javapipeline.PriceClusterReads;
import gedi.riboseq.javapipeline.PriceCodonInference;
import gedi.riboseq.javapipeline.PriceCodonViewerIndices;
import gedi.riboseq.javapipeline.PriceCollectSufficientStatistics;
import gedi.riboseq.javapipeline.PriceDetermineDelta;
import gedi.riboseq.javapipeline.PriceEstimateModel;
import gedi.riboseq.javapipeline.PriceIdentifyMaxPos;
import gedi.riboseq.javapipeline.PriceMultipleTestingCorrection;
import gedi.riboseq.javapipeline.PriceNoiseTraining;
import gedi.riboseq.javapipeline.PriceOptimisticCodonMapping;
import gedi.riboseq.javapipeline.PriceOrfInference;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.riboseq.javapipeline.PriceReassignCodons;
import gedi.riboseq.javapipeline.PriceSetupOrfInference;
import gedi.riboseq.javapipeline.PriceSignalToNoise;
import gedi.riboseq.javapipeline.PriceStartPredictionTraining;
import gedi.riboseq.javapipeline.PriceWriteCodons;
import gedi.util.LogUtils.LogMode;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class Price {

	
	public static void main(String[] args) {
		System.setProperty("smile.threads", "1");
		
		
		Gedi.startup(true,LogMode.Normal,"Price");
		
		PriceParameterSet params = new PriceParameterSet();
		GediProgram pipeline = GediProgram.create("PRICE",
				new PriceIdentifyMaxPos(params),
				new PriceCollectSufficientStatistics(params),
				new PriceEstimateModel(params),
				new PriceClusterReads(params),
				new PriceSetupOrfInference(params),
				new PriceCodonInference(params),
				new PriceDetermineDelta(params),
				new PriceWriteCodons(params),
				new PriceCodonViewerIndices(params),
				new PriceStartPredictionTraining(params),
				new PriceNoiseTraining(params),
				new PriceOrfInference(params),
				new PriceMultipleTestingCorrection(params),
				new PriceReassignCodons(params),
				new PriceSignalToNoise(params),
				new PriceOptimisticCodonMapping(params)
				);
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("PRICE","PRICE is an analysis method for Ribo-seq data.",args));
		
	}
	
}
