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
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.riboseq.inference.orf.PriceOrfType;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.startup.PriceStartUp;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.SparseMemoryFloatArray;
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

public class UOrfReadsPriceAnalysis implements PriceAnalysis<Void> {

	/**
	 * For auto inclusion in {@link PriceStartUp}
	 */
	public static final String name = "uorfreads";
	public static final Function<PriceParameterSet,GediParameter[]> params = set->new GediParameter[] {
			set.indices
	};
	
	
	private CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codons;
	
	public UOrfReadsPriceAnalysis(String[] conditions, PriceParameterSet param) throws IOException {
		codons = new CenteredDiskIntervalTreeStorage<>(param.indices.getFile().getPath());
	}
	
	@Override
	public void process(MajorIsoform data, LineWriter out, Void ctx) {
		
		GenomicRegion utrreg = data.getTranscript().map(new ArrayGenomicRegion(0,data.getOrf(0).getStart()));
		double total = codons.ei(new ImmutableReferenceGenomicRegion<>(data.getTranscript().getReference(), utrreg))
			.filter(c->utrreg.containsUnspliced(c.getRegion()))
			.mapToDouble(c->c.getData().sum()).sum();
		
		out.writef2("%s\t%.1f\n", data.getTranscript().getData().getGeneId(),
				total
				);
	}

	@Override
	public void header(LineWriter out) throws IOException {
		out.writeLine("GeneId\tuORF reads");
	}

	@Override
	public Void createContext() {
		return null;
	}

	@Override
	public void reduce(List<Void> ctx, LineWriter out) {
	}

	
	
}
