package gedi.macoco.javapipeline;

import java.io.File;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.StringUtils;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.inference.isoforms.EquivalenceClassEffectiveLengths;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.userInteraction.progress.Progress;

public class EquivalenceClassesEffectiveLengths extends GediProgram {



	public EquivalenceClassesEffectiveLengths(MacocoParameterSet params) {
		
		addInput(params.countTable);
		addInput(params.lenDistTable);
		addInput(params.genomic);
		addInput(params.mrnas);
		addInput(params.strandness);
		
		addInput(params.prefix);
		addOutput(params.elTable);
	}

	

	@Override
	public String execute(GediProgramContext context) throws Exception {
		File countsf = getParameter(0);
		File lenf= getParameter(1);
		Genomic genomic = getParameter(2);
		GenomicRegionStorage<NameProvider> mRNAs = getParameter(3);
		Strandness strand = getParameter(4);
		
		context.getLog().info("Running "+StringUtils.removeFooter(StringUtils.removeHeader(getOutputSpec().get(0).getName(),"${prefix}."),".tsv")+"...");
		
		
		MemoryIntervalTreeStorage<Transcript> trans;
		if (mRNAs==null) {
			trans = genomic.getTranscripts();
		} else {
			trans = new MemoryIntervalTreeStorage<>(Transcript.class);
			trans.fill(mRNAs.ei().map(r->new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(),new Transcript(r.getData().getName(), r.getData().getName(), -1, -1))));
		}

		EquivalenceClassInfo info = new EquivalenceClassInfo(countsf.getAbsolutePath(), lenf.getAbsolutePath(), null);
		
		double[][] el = new double[info.getConditionNames().length][info.getE().length];
		for (int cond=0; cond<el.length; cond++) {
			double[] eff = EquivalenceClassEffectiveLengths.preprocessEff(info.getReadLengths()[cond]);
			
			
			ExtendedIterator<ImmutableReferenceGenomicRegion<String>> tit = trans.ei().progress(context.getProgress(), (int)trans.size(), (x)->"Collecting equivalence classes")
					.map(r->{
						ReferenceSequence ref = r.getReference();
						if (strand==Strandness.Unspecific)
							ref = ref.toStrandIndependent();
						else if (strand==Strandness.Antisense)
							ref = ref.toOppositeStrand();
						return new ImmutableReferenceGenomicRegion<String>(ref, r.getRegion(), r.getData().getTranscriptId());
					});
			EquivalenceClassEffectiveLengths<String> len = new EquivalenceClassEffectiveLengths<String>(tit,eff);
			
			Progress prog = context.getProgress();
			prog.init().setCount(info.getE().length).setDescription("Processing "+(cond+1)+"/"+el.length);
			for (int i=0; i<info.getE().length; i++) {
				el[cond][i] = len.getEffectiveLength(info.getE()[i]);
				prog.incrementProgress();
			}
			prog.finish();
		}
		
		
		LineWriter out = new LineOrientedFile(getOutputFile(0).getPath()).write();
		out.writef("Equivalence class");
		for (int i=0; i<info.getConditionNames().length; i++)
			out.writef("\t%s", info.getConditionNames()[i]);
		out.writeLine();
		for (int i=0; i<info.getE().length; i++) {
			out.writef(StringUtils.concat(",", info.getE()[i]));
			for (int cond=0; cond<el.length; cond++)
				out.writef("\t%.2f", el[cond][i]);
			out.writeLine();
		}
		out.close();
		
		return null;
	}


}
