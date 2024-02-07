package gedi.riboseq.javapipeline;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Predicate;
import java.util.logging.Level;

import gedi.app.extension.ExtensionContext;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.numeric.diskrmq.DiskGenomicNumericBuilder;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStorageCapabilities;
import gedi.core.region.GenomicRegionStorageExtensionPoint;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.sequence.SequenceProvider;
import gedi.riboseq.cleavage.RiboModel;
import gedi.riboseq.inference.clustering.RiboClusterInfo;
import gedi.riboseq.inference.codon.Codon;
import gedi.riboseq.inference.codon.CodonInference;
import gedi.riboseq.inference.orf.NoiseModel;
import gedi.riboseq.inference.orf.OrfInference;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.riboseq.inference.orf.StartCodonScorePredictor;
import gedi.riboseq.inference.orf.StartCodonTraining;
import gedi.riboseq.utils.RiboUtils;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.MemoryFloatArray;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.tree.redblacktree.IntervalTreeSet;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.testing.MultipleTestingCorrection;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;

public class PriceMultipleTestingCorrection extends GediProgram {

	public PriceMultipleTestingCorrection(PriceParameterSet params) {
		addInput(params.prefix);
		addInput(params.pvals);
		addInput(params.orfinference);
		addInput(params.fdr);
		addInput(params.orfsbin);
		addInput(params.keepAnno);
		
		addOutput(params.orfstmp);
	}
	
	public String execute(GediProgramContext context) throws IOException {
		
		String prefix = getParameter(0);
		File pvalFile = getParameter(1);
		OrfInference v = getParameter(2);
		double fdr = getParameter(3);
		boolean keepAnno = getBooleanParameter(5);
		
		DoubleArrayList pvals = new DoubleArrayList();
		PageFile pvf = new PageFile(pvalFile.getAbsolutePath());
		while (!pvf.eof())
			pvals.add(pvf.getDouble());
		pvf.close();
		
		
		context.getLog().log(Level.INFO, "Found "+pvals.size()+" ORFs");
		context.getLog().log(Level.INFO, "Multiple testing correction and filtering");
		double[] corr = MultipleTestingCorrection.benjaminiHochberg(pvals.toDoubleArray());
		PageFile in = new PageFile(prefix+".orfs.bin");
		in.getContext().add(Class.class, PriceOrf.class);
		MutableInteger index = new MutableInteger();
		MutableInteger fdrn = new MutableInteger(0);
		
		GenomicRegionStorage out = GenomicRegionStorageExtensionPoint.getInstance().get(new ExtensionContext().add(String.class, prefix+".tmp.orfs").add(Class.class, PriceOrf.class), GenomicRegionStorageCapabilities.Disk, GenomicRegionStorageCapabilities.Fill);
		out.fill(
				in.ei().map(pf->{
					try {
						MutableReferenceGenomicRegion<PriceOrf> orf = new MutableReferenceGenomicRegion<>();
						orf.deserialize(in);
						v.setCorrectedPvalue(orf,corr[index.N++]);
						return orf;
					} catch (IOException e) {
						throw new RuntimeException("Could not read temporary orfs!",e);
					}
					})
				.progress(context.getProgress(), corr.length, r->r.toLocationStringRemovedIntrons())
				.filter(o->o.getData().getCombinedP()<fdr || (keepAnno && v.isAnnotatedEnd(o, true)))
				.sideEffect(o->{
					if (o.getData().getCombinedP()<fdr) fdrn.N++;
				})
				.iff(v.getCheckOrfNames().size()>0, ei->ei.sideEffect(o->v.setDetected(o)))
				);
		in.close();
		
		context.getLog().log(Level.INFO, "Remaining after multiple testing correction: "+fdrn.N+" ORFs");
		if (keepAnno)
			context.getLog().log(Level.INFO, "+ remaining annotated: "+(out.size()-fdrn.N)+" ORFs");
		
		
		return null;
	}
	

}
