package gedi.riboseq.analysis;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.startup.PriceStartUp;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.collections.PositionIterator;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.distributions.LfcDistribution;
import gedi.util.mutable.MutablePair;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;
import gedi.util.r.RRunner;

public class UtrLengthPriceAnalysis implements PriceAnalysis<Void> {

	/**
	 * For auto inclusion in {@link PriceStartUp}
	 */
	public static final String name = "transcript";

	
	
	@Override
	public void process(MajorIsoform data, LineWriter out, Void ctx) {
		out.writef2("%s\t%s\n", data.getTranscript().getData().getGeneId(),data.getTranscript().getData().getTranscriptId());
	}

	@Override
	public void header(LineWriter out) throws IOException {
		out.writeLine("GeneId\tTranscriptId");
	}

	@Override
	public Void createContext() {
		return null;
	}

	@Override
	public void reduce(List<Void> ctx, LineWriter out) {
	}

	
	
}
