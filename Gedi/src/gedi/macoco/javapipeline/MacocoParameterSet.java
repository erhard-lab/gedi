package gedi.macoco.javapipeline;

import java.io.File;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.feature.special.Downsampling;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.LongParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class MacocoParameterSet extends GediParameterSet {

	
	public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
	
	public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
	public GediParameter<Long> seed = new GediParameter<Long>(this,"seed", "Random number generator seed", false, new LongParameterType(), 42L);
	
	public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "The mapped reads from the ribo-seq experiment.", false, new StorageParameterType<AlignedReadsData>());
	public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"genomic", "The indexed GEDI genome.", true, new GenomicParameterType());
	public GediParameter<GenomicRegionStorage<NameProvider>> mrnas = new GediParameter<GenomicRegionStorage<NameProvider>>(this,"mRNA", "Use other mRNAs than specified by the genomic object.", false, new StorageParameterType<NameProvider>(),true);
	public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Which strandness.", false, new EnumParameterType<>(Strandness.class));
	public GediParameter<Downsampling> down = new GediParameter<Downsampling>(this,"downsampling", "Which downsampling.", false, new EnumParameterType<>(Downsampling.class),Downsampling.No);
	
	public GediParameter<File> countTable = new GediParameter<File>(this,"${prefix}.counts.tsv", "File containing the counts table", false, new FileParameterType());
	public GediParameter<File> lenDistTable = new GediParameter<File>(this,"${prefix}.lengths.tsv", "File containing the length distribution", false, new FileParameterType());
	public GediParameter<File> outTable = new GediParameter<File>(this,"${prefix}.tsv", "File containing the output table", false, new FileParameterType());
	public GediParameter<File> emMLTable = new GediParameter<File>(this,"${prefix}.ML.tsv", "File containing the tpms estimated by the EM", false, new FileParameterType());
	public GediParameter<File> lassoTable = new GediParameter<File>(this,"${prefix}.lasso.tsv", "File containing the tpms estimated by the penalized EM", false, new FileParameterType());
	public GediParameter<File> macocoTable = new GediParameter<File>(this,"${prefix}.macoco.tsv", "File containing the tpms estimated by the MACOCO", false, new FileParameterType());
	public GediParameter<File> elTable = new GediParameter<File>(this,"${prefix}.effLen.tsv", "File containing the effective lengths of all equivalence classes", false, new FileParameterType());
	
	
	
}
