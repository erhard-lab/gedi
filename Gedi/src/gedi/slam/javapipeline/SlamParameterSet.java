package gedi.slam.javapipeline;

import java.io.File;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.ChoicesParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StorageParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class SlamParameterSet extends GediParameterSet {

	public GediParameter<File> paramFile = new GediParameter<File>(this,"${prefix}.param", "File containing the parameters used to call PRICE", false, new FileParameterType());
	public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());
	
	public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "The prefix used for all output files", false, new StringParameterType());
	public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "A test", false, new BooleanParameterType());
	public GediParameter<File> testFile = new GediParameter<File>(this,"${prefix}.test", "Output file for test", false, new FileParameterType());
	
	public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
	public GediParameter<GenomicRegionStorage<AlignedReadsData>> reads = new GediParameter<GenomicRegionStorage<AlignedReadsData>>(this,"reads", "The mapped reads from the SLAM-seq experiment.", false, new StorageParameterType<AlignedReadsData>());
	public GediParameter<GenomicRegionStorage<NameProvider>> locations = new GediParameter<GenomicRegionStorage<NameProvider>>(this,"locations", "Map reads to special locations in the genome.", false, new StorageParameterType<NameProvider>(),true);
	public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"genomic", "The indexed GEDI genome.", true, new GenomicParameterType());
	public GediParameter<Boolean> full = new GediParameter<Boolean>(this,"full", "Full output (including credible interval, coverages etc.", false, new BooleanParameterType());
	public GediParameter<Boolean> newsnp = new GediParameter<Boolean>(this,"newsnp", "New way of estimating snps", false, new BooleanParameterType());
	
	public GediParameter<Integer> trim5p = new GediParameter<Integer>(this,"trim5p", "The number bases to trim at the 5' ends of reads", false, new IntParameterType(), 0);
	public GediParameter<Integer> trim3p = new GediParameter<Integer>(this,"trim3p", "The number bases to trim at the 3' ends of reads", false, new IntParameterType(), 0);
	public GediParameter<Strandness> strandness = new GediParameter<Strandness>(this,"strandness", "Whether sequencing protocol was stranded (Sensse), strand unspecific (Unspecific), or opposite strand (Antisense).", false, new EnumParameterType<>(Strandness.class),true);
	public GediParameter<String> viral = new GediParameter<String>(this,"viral", "Name of the viral chromosome to count mismatches from.", false, new StringParameterType(),true);
	public GediParameter<Boolean> all = new GediParameter<Boolean>(this,"allGenes", "Keep all genes", false, new BooleanParameterType());
	public GediParameter<Boolean> sparseOutput = new GediParameter<Boolean>(this,"sparse", "write sparse matrix format", false, new BooleanParameterType());
	public GediParameter<Boolean> m10x = new GediParameter<Boolean>(this,"10x", "equivalent to -sparse -lenient -strandness Sense -conv new", false, new BooleanParameterType()).setShortcut("-sparse","-lenient","-strandness","Sense","-conv","new");
	public GediParameter<Boolean> doubleOnly = new GediParameter<Boolean>(this,"double", "Use only doubly sequenced parts for estimation", false, new BooleanParameterType());
	public GediParameter<Boolean> modelall = new GediParameter<Boolean>(this,"modelall", "Use all reads (inside of a gene) for estimating the models", false, new BooleanParameterType());
	public GediParameter<Boolean> highmem = new GediParameter<Boolean>(this,"highmem", "Do not stream reads, but load all for a gene into memory", false, new BooleanParameterType());
	

	public GediParameter<File> snpFile = new GediParameter<File>(this,"${prefix}.snpdata", "File containing the SNP data", false, new FileParameterType()).setRemoveFile(true);
	public GediParameter<File> strandnessFile = new GediParameter<File>(this,"${prefix}.strandness", "File containing strandness information", false, new FileParameterType());
	
	public GediParameter<File> dataFile = new GediParameter<File>(this,"${prefix}.genedata", "File containing the gene data", false, new FileParameterType()).setRemoveFile(true);
	public GediParameter<File> dataExtFile = new GediParameter<File>(this,"${prefix}.extdata", "File containing the summary data", false, new FileParameterType()).setRemoveFile(true);
	public GediParameter<File> mismatchFile = new GediParameter<File>(this,"${prefix}.mismatches.tsv", "File containing the mismatch rates", false, new FileParameterType());
	public GediParameter<File> mismatchPosFile = new GediParameter<File>(this,"${prefix}.mismatchdetails.tsv", "File containing the mismatch details (positions, overlap, ...)", false, new FileParameterType());
	public GediParameter<File> binomFile = new GediParameter<File>(this,"${prefix}.binom.tsv", "File containing the binomial mixture statistics", false, new FileParameterType());
	public GediParameter<File> binomOverlapFile = new GediParameter<File>(this,"${prefix}.binomOverlap.tsv", "File containing the binomial mixture statistics for overlaps", false, new FileParameterType());
	public GediParameter<File> ntrFile = new GediParameter<File>(this,"${prefix}.ntrstat.tsv", "File containing NTR statistics", false, new FileParameterType());
	
	public GediParameter<Boolean> numi = new GediParameter<Boolean>(this,"nUMI", "Compute nUMIs", false, new BooleanParameterType());
	public GediParameter<Boolean> introns = new GediParameter<Boolean>(this,"introns", "Compute introns", false, new BooleanParameterType());
	public GediParameter<Boolean> lenientOverlap = new GediParameter<Boolean>(this,"lenient", "Apply lenient transcript compatibility (the two leading and trailing nt may be intronic, reads may extend beyond transcript boundaries); use this for 10x CellRanger mapped reads!", false, new BooleanParameterType());
	public GediParameter<Double> biasnew = new GediParameter<Double>(this,"bias", "Also produce output for biased estimates of p_conv (by this factor)", false, new DoubleParameterType(),0.0);
	
//	public GediParameter<File> overlapFile = new GediParameter<File>(this,"${prefix}.overlap.tsv", "File containing the overlap conversion statistics", false, new FileParameterType());
	public GediParameter<File> doublehitFile = new GediParameter<File>(this,"${prefix}.doublehit.tsv", "File containing the double hit statistics", false, new FileParameterType());
	public GediParameter<File> binomOut = new GediParameter<File>(this,"${prefix}.binomEstimated.tsv", "File with estimation statistics", false, new FileParameterType());
	public GediParameter<File> binomOverlapOut = new GediParameter<File>(this,"${prefix}.binomOverlapEstimated.tsv", "File with estimation statistics", false, new FileParameterType());
	public GediParameter<File> errorOut = new GediParameter<File>(this,"${prefix}.errorEstimated.tsv", "File with estimation statistics", false, new FileParameterType());
	
	public GediParameter<Integer> minEstimateReads = new GediParameter<Integer>(this,"minEstimateReads", "Specify the minimal number of reads to be used for parameter estimation!", false, new IntParameterType(),10000);

	public GediParameter<File> out10x = new GediParameter<File>(this,"${prefix}.matrix", "Directory containing the output matrix", false, new FileParameterType());
	public GediParameter<File> outTable = new GediParameter<File>(this,"${prefix}.tsv.gz", "File containing the output table", false, new FileParameterType());
	public GediParameter<File> extTable = new GediParameter<File>(this,"${prefix}.ext.tsv", "File containing the summary output table", false, new FileParameterType());
	public GediParameter<File> rateTable = new GediParameter<File>(this,"${prefix}.rates.tsv", "File containing the output table", false, new FileParameterType());
	
	public GediParameter<String> conv = new GediParameter<String>(this,"conv", "Specific the conversion rate", false, new StringParameterType(),true);
	public GediParameter<String> err = new GediParameter<String>(this,"err", "Specific the error rate instead of estimate it from no4sU", false, new StringParameterType(),true);
	public GediParameter<String> errlm = new GediParameter<String>(this,"errlm", "Specify the error rate correction linear model (which is applied to the other 11 error rates; expression is evaluated by JS, variables are AC, AG, AT, ... and median)", false, new StringParameterType(),"TA*1.434");//,"TA*0.68245+TG*0.00796");
	public GediParameter<String> errlm2 = new GediParameter<String>(this,"errlm2", "Specify the error rate correction linear model (which is applied to the other 11 error rates; expression is evaluated by JS, variables are AC, AG, AT, ... and median)", false, new StringParameterType(),"TA*1.434");//,"TA*0.68245+TG*0.00796");
	public GediParameter<String> no4sUpattern = new GediParameter<String>(this,"no4sUpattern", "The pattern to recognize no4sU conditions", false, new StringParameterType(),"no4sU|nos4U");
	
	public GediParameter<Double> lower = new GediParameter<Double>(this,"lower", "The lower credible interval bound", false, new DoubleParameterType(), 0.05);
	public GediParameter<Double> upper = new GediParameter<Double>(this,"upper", "The upper credible interval bound", false, new DoubleParameterType(), 0.95);
	
	public GediParameter<Double> snpConv = new GediParameter<Double>(this,"snpConv", "Conversion rate for SNP calling", false, new DoubleParameterType(), 0.3);
	public GediParameter<Double> snpPval = new GediParameter<Double>(this,"snppval", "Minimal posterior probability of being a conversion", false, new DoubleParameterType(), 0.001);
	
	
	public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Use R to produce various plots",false, new BooleanParameterType());
	
	public GediParameter<ReadCountMode> mode = new GediParameter<ReadCountMode>(this,"mode", "Read count mode what to do with multi-mapping reads", false, new ChoicesParameterType<ReadCountMode>(ReadCountMode.valueOf, ReadCountMode.class),ReadCountMode.Weight);
	public GediParameter<ReadCountMode> overlap = new GediParameter<ReadCountMode>(this,"overlap", "Overlapping gene mode: What to do for locations that are compatible with more than one gene", false, new ChoicesParameterType<ReadCountMode>(ReadCountMode.valueOf, ReadCountMode.class),ReadCountMode.All);

	
	
}
