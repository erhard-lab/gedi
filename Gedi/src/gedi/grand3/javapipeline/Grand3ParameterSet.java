package gedi.grand3.javapipeline;

import java.io.File;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.grand3.estimation.estimators.ModelEstimationMethod;
import gedi.grand3.javapipeline.Grand3SetupTargetsGenes.ExonIntronMode;
import gedi.grand3.processing.Resimulator.ResimulatorModelType;
import gedi.grand3.targets.TargetCollection;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.ChoicesParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.InternalParameterType;
import gedi.util.program.parametertypes.LongParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class Grand3ParameterSet extends GediParameterSet {

	public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call Grand3", false, new FileParameterType());
	public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());
	
	public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
	
	public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
	public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "The mapped reads from the SLAM-seq experiment.", false, new StorageParameterType<AlignedReadsData>());
	public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"genomic", "The indexed GEDI genome.", true, new GenomicParameterType());
	public GediParameter<Boolean> noplot = new GediParameter<Boolean>(this,"noplot", "Don't use R to produce various during runtime (they can be produced afterwars by running the generated R script)",false, new BooleanParameterType());
	
	public GediParameter<Boolean> nosnp = new GediParameter<Boolean>(this,"nosnps", "Do not consider any snps", false, new BooleanParameterType());
	public GediParameter<String> snp = new GediParameter<String>(this,"snps", "File containing precomputed SNPs", false, new StringParameterType(),true);
	public GediParameter<Double> snpConv = new GediParameter<Double>(this,"snpConv", "Conversion rate for SNP calling", false, new DoubleParameterType(), 0.3);
	public GediParameter<Double> snpPval = new GediParameter<Double>(this,"snppval", "Minimal posterior probability of being a conversion", false, new DoubleParameterType(), 0.001);
	public GediParameter<Boolean> blacklistSnp = new GediParameter<Boolean>(this,"blacklistControlMM", "Blacklist all mismatches observed in the (no4sU) control", false, new BooleanParameterType());
	
	public GediParameter<String> clip = new GediParameter<String>(this,"clip", "Force clipping (allowed: none/auto/5pr1,5pr2)", false, new StringParameterType(), "auto");
	public GediParameter<String> subreadsToUse = new GediParameter<String>(this,"restrict-subreads", "Only use given subreads for estimation (of genewise parameters, specify indices or names, separated by space)", true, new StringParameterType(),true);
	public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Whether sequencing protocol was stranded (Sensse), strand unspecific (Unspecific), or opposite strand (Antisense).", false, new EnumParameterType<>(Strandness.class),true);
	public GediParameter<File> strandnessFile = new GediParameter<File>(this,"${prefix}.strandness.tsv", "Inferred strandness file", false, new FileParameterType());

	public GediParameter<ReadCountMode> mode = new GediParameter<ReadCountMode>(this,"mode", "Read count mode what to do with multi-mapping reads", false, new ChoicesParameterType<ReadCountMode>(ReadCountMode.valueOf, ReadCountMode.class),ReadCountMode.Weight);
	public GediParameter<ReadCountMode> overlap = new GediParameter<ReadCountMode>(this,"overlap", "Overlapping gene mode: What to do for locations that are compatible with more than one gene", false, new ChoicesParameterType<ReadCountMode>(ReadCountMode.valueOf, ReadCountMode.class),ReadCountMode.Unique);
	public GediParameter<Integer> introntol = new GediParameter<Integer>(this,"introntol", "Tolerance for introns", false, new IntParameterType(), 5);
	
	public GediParameter<TargetCollection> targetCollection = new GediParameter<TargetCollection>(this,"target-description","Internal use only",false,new InternalParameterType<>(TargetCollection.class));
	public GediParameter<String> targets = new GediParameter<String>(this,"targets","Define which targets to estimate parameters for",false,new StringParameterType(),"genes");
	public GediParameter<ExonIntronMode> genemode = new GediParameter<ExonIntronMode>(this, "genemode", "Specify whether exons/intron should be used for (i) estimating global parameters and (ii) estimating NTRs", false, new EnumParameterType<>(ExonIntronMode.class));

	public GediParameter<Boolean> debug = new GediParameter<Boolean>(this,"debug", "Debug mode (don't use it!)",false, new BooleanParameterType()).setShortcut("-nthreads","0");
	public GediParameter<Boolean> subreadCit = new GediParameter<Boolean>(this,"subread", "Write subread cit file",false, new BooleanParameterType());
	public GediParameter<Boolean> burstMCMC = new GediParameter<Boolean>(this,"burstMCMC", "Write burstMCMC input file",false, new BooleanParameterType());
	
	public GediParameter<Boolean> forceDense = new GediParameter<Boolean>(this,"output-table", "Force output as a (dense) table (default if # conditions <=30)",false, new BooleanParameterType());
	public GediParameter<Boolean> forceSparse = new GediParameter<Boolean>(this,"output-sparse", "Force output as a sparse matrix (default if # conditions >30)",false, new BooleanParameterType());
	
	
	public GediParameter<File> snpFile = new GediParameter<File>(this,"${prefix}.snpdata.gz", "File containing the SNP data", false, new FileParameterType());
	public GediParameter<File> mismatchPositionFile = new GediParameter<File>(this,"${prefix}.mismatch.raw.position.tsv.gz", "How many mismatches remain after correction due to overlapping reads", false, new FileParameterType());
	
	public GediParameter<Integer> clip5p1 = new GediParameter<Integer>(this,"clip5p1", "The number bases to clip at the 5' ends of read 1", false, new IntParameterType(), -1);
	public GediParameter<Integer> clip3p1 = new GediParameter<Integer>(this,"clip3p1", "The number bases to clip at the 3' ends of read 1", false, new IntParameterType(), -1);
	public GediParameter<Integer> clip5p2 = new GediParameter<Integer>(this,"clip5p2", "The number bases to clip at the 5' ends of read 2", false, new IntParameterType(), -1);
	public GediParameter<Integer> clip3p2 = new GediParameter<Integer>(this,"clip3p2", "The number bases to clip at the 3' ends of read 2", false, new IntParameterType(), -1);
	public GediParameter<File> clipFile = new GediParameter<File>(this,"${prefix}.clip.tsv", "Inferred clipping parameters", false, new FileParameterType());

	
	public GediParameter<File> mismatchFile = new GediParameter<File>(this,"${prefix}.conversion.freq.tsv.gz", "File containing the mismatch frequencies", false, new FileParameterType());
	public GediParameter<File> perrFile = new GediParameter<File>(this,"${prefix}.conversion.perr.tsv.gz", "File containing the mismatch frequencies for perr estimation", false, new FileParameterType());
	public GediParameter<File> knMatFile = new GediParameter<File>(this,"${prefix}.conversion.knmatrix.tsv.gz", "File containing the sufficient statistics for estimating the mixture model", false, new FileParameterType());
	public GediParameter<File> subreadSemanticFile = new GediParameter<File>(this,"${prefix}.subread.tsv", "File containing the semantics of subreads", false, new FileParameterType());
	public GediParameter<File> experimentalDesignFile = new GediParameter<File>(this,"${prefix}.experimentalDesign.tsv", "File containing the label design", false, new FileParameterType());
	public GediParameter<File> allPseudobulkFile = new GediParameter<File>(this,"${prefix}.pseudobulk.all.tsv.gz", "Pseudobulk definition file with all cells from the same sample being a pseudobulk; can be used to call grand3 again and using this as -pseudobulkFile parameter", false, new FileParameterType());
	
	public GediParameter<Boolean> auto = new GediParameter<Boolean>(this,"auto", "Do not stop analysis if there is no experimental design (i.e. try to infer everything from the sample names!",false, new BooleanParameterType());
	public GediParameter<Boolean> confidence = new GediParameter<Boolean>(this,"ci", "Compute 95% confidence intervals via the profile likelihood!",false, new BooleanParameterType());
	public GediParameter<Boolean> profile = new GediParameter<Boolean>(this,"profile", "Compute and output the likelihood profiles (implies -ci)!",false, new BooleanParameterType());
	
	public GediParameter<Boolean> resim = new GediParameter<Boolean>(this,"resim", "Resimulate reads",false, new BooleanParameterType());
	public GediParameter<ResimulatorModelType> resimModelType = new GediParameter<ResimulatorModelType>(this,"resim-model", "Which model for resimulate (default: Do not resimulate, but write the reads as is).",false, new EnumParameterType<ResimulatorModelType>(ResimulatorModelType.class),true);
	public GediParameter<Long> resimSeed = new GediParameter<Long>(this,"resim-seed", "Seed for resimulation.",false, new LongParameterType(),1337L);
	public GediParameter<File> resimFile = new GediParameter<File>(this,"${prefix}.resim.${resim-model}.cit", "File containing the resimulated reads", false, new FileParameterType());
	
	public GediParameter<File> burstMCMCFile = new GediParameter<File>(this,"${prefix}.burstMCMC.txt", "File containing the burstMCMC inputs", false, new FileParameterType());
	
	public GediParameter<File> modelFile = new GediParameter<File>(this,"${prefix}.model.parameters.tsv", "File containing the mixture model parameters", false, new FileParameterType());
	public GediParameter<File> modelBinFile = new GediParameter<File>(this,"${prefix}.model.parameters.bin", "File containing the mixture model parameters", false, new FileParameterType());
	
	public GediParameter<ModelEstimationMethod> estimMethod = new GediParameter<ModelEstimationMethod>(this, "estim", "Which method to use for estimating the global parameters. If more than one are given, the first is used to estimate target parameters", true, new EnumParameterType<>(ModelEstimationMethod.class),ModelEstimationMethod.Full);
	
	public GediParameter<String> targetMixmat = new GediParameter<String>(this,"targetMixmat", "Which target to output the MixMatrix for", false, new StringParameterType(),true);
	public GediParameter<File> targetMixmatFile = new GediParameter<File>(this,"${prefix}.mixmatrices/${targetMixmat}.knmatrix.tsv.gz", "Mixmatrix for the given target", false, new FileParameterType());
	
	public GediParameter<String> targetMergeTable = new GediParameter<String>(this,"targetMergeTab", "A tsv file containing columns name (original name) and merged (new name)", false, new StringParameterType(),true);
	
	public GediParameter<String> targetsName = new GediParameter<String>(this,"targetsName", "Name for the target output files", false, new StringParameterType(),"targets");
	public GediParameter<File> targetBinFile = new GediParameter<File>(this,"${prefix}.${targetsName}.bin", "File containing the raw estimates per target", false, new FileParameterType(true));
	public GediParameter<File> targetFolder = new GediParameter<File>(this,"${prefix}.${targetsName}", "Folder containing sparse matrices of all estimates per target", false, new FileParameterType());

	public GediParameter<String> pseudobulkFile = new GediParameter<String>(this,"pseudobulkFile", "File containing merge and distribute information for output (target) conditions (cells), e.g. useful for creating pseudobulk samples; if specified, no normal target output is generated.", false, new StringParameterType(),true);
	public GediParameter<String> pseudobulkName = new GediParameter<String>(this,"pseudobulkName", "The prefix used for pseudobulk output", false, new StringParameterType(),"all");
	public GediParameter<Double> pseudobulkMinimalPurity = new GediParameter<Double>(this,"pseudobulkPurity", "The minimal purity (how many cells are from a single library/sample) of a pseudobulk sample to estimate the parameters.", false, new DoubleParameterType(),0.95);
	
	public GediParameter<File> pseudobulkBinFile = new GediParameter<File>(this,"${prefix}.pseudobulk.${targetsName}.${pseudobulkName}.bin", "File containing the raw estimates per pseudobulk target", false, new FileParameterType(true)) ;
	public GediParameter<File> pseudobulkFolder = new GediParameter<File>(this,"${prefix}.pseudobulk.${targetsName}.${pseudobulkName}", "Folder containing sparse matrices of all estimates per pseudobulk target", false, new FileParameterType());

	public GediParameter<String> jointModel = new GediParameter<String>(this,"jointModel", "Definition of which models to estimate jointly (done on the . separated name, i.e. this should be a 1010110 with as many digits as . separated fields, and samples with the same remaining name after removing 0 fields are joined)", false, new StringParameterType(),true);

	
	public GediParameter<String> test = new GediParameter<String>(this,"test", "Which gene or location for testing?", false, new StringParameterType(),true).setShortcut("-nthreads","0");
		
	
	

	
	
	
}
