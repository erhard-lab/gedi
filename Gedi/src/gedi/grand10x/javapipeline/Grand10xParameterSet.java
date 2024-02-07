package gedi.grand10x.javapipeline;

import java.io.File;

import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionStorage;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class Grand10xParameterSet extends GediParameterSet {

	public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
	public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "Output prefix", true, new StringParameterType());
	public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"g", "Genomic name", true, new GenomicParameterType());
	public GediParameter<GenomicRegionStorage<BarcodedAlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<BarcodedAlignedReadsData>>(this,"reads", "The mapped reads from the 10x experiment.", false, new StorageParameterType<BarcodedAlignedReadsData>());
	public GediParameter<Integer> maxdist = new GediParameter<Integer>(this,"maxdist", "The maximal distance between two read 5' ends", false, new IntParameterType(), 200);
	public GediParameter<String> whitelist = new GediParameter<String>(this,"cells", "File containing the cell barcodes to use", false, new StringParameterType());
	public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Produce plots", false, new BooleanParameterType());

	public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "Only search against chromosome 22", false, new BooleanParameterType());
	public GediParameter<Boolean> computeStat = new GediParameter<Boolean>(this,"stat", "Compute mismatch statistics", false, new BooleanParameterType());
	public GediParameter<Integer> mincov = new GediParameter<Integer>(this,"mincov", "Minimum coverage for output", false, new IntParameterType(), 0);
	
	
	public GediParameter<File> clustercit = new GediParameter<File>(this,"${prefix}.cluster.cit", "Clusters for the viewer", false, new FileParameterType());
	public GediParameter<File> umicit = new GediParameter<File>(this,"${prefix}.umi.cit", "Umi reads", false, new FileParameterType());
	public GediParameter<File> infer3ptsv = new GediParameter<File>(this,"${prefix}.3p.tsv", "Table of 3' end information.", false, new FileParameterType());
	public GediParameter<File> infer3pparam = new GediParameter<File>(this,"${prefix}.3p.parameter", "3' parameter file.", false, new FileParameterType());
	
	public GediParameter<Double> snpConv = new GediParameter<Double>(this,"snpConv", "Conversion rate for SNP calling", false, new DoubleParameterType(), 0.2);
	public GediParameter<Double> snpPval = new GediParameter<Double>(this,"snppval", "Minimal posterior probability of being a conversion", false, new DoubleParameterType(), 0.01);

}