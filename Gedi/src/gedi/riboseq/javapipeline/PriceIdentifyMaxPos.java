package gedi.riboseq.javapipeline;


import java.io.File;
import java.util.function.Supplier;
import java.util.logging.Logger;

import javax.script.ScriptException;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.feature.AlignedReadsDataToFeatureProgram;
import gedi.core.region.feature.GenomicRegionFeatureProgram;
import gedi.core.region.feature.features.AbsolutePosition;
import gedi.core.region.feature.features.AnnotationFeature;
import gedi.core.region.feature.output.FeatureStatisticOutput;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.dataframe.DataFrame;
import gedi.util.io.text.tsv.formats.Csv;
import gedi.util.program.GediParameter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.userInteraction.progress.Progress;

public class PriceIdentifyMaxPos extends GediProgram {

	
	
	public PriceIdentifyMaxPos(PriceParameterSet params) {
		addInput(params.prefix);
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.reads);
		addOutput(params.maxPos);
	}
	
	
	public String execute(GediProgramContext context) throws ScriptException {
		
		String prefix= getParameter(0);
		int nthreads = getIntParameter(1);
		Genomic genomic = getParameter(2);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(3);
		
		context.getLog().info("Estimating maxpos...");
		GenomicRegionFeatureProgram<AlignedReadsData> program = new GenomicRegionFeatureProgram<AlignedReadsData>();
		program.setThreads(nthreads);
		
		AnnotationFeature<Transcript> a = new AnnotationFeature<Transcript>(false);
		a.addTranscripts(genomic);
		a.setId("transcript");
		program.add(a);
		
		AbsolutePosition p = new AbsolutePosition();
		p.setReportFurtherDownstream(false);
		p.setReportFurtherUpstream(false);
		p.setAnnotationPosition(GenomicRegionPosition.StartCodon);
		p.setId("Position");
		program.add(p,"transcript");
		
		FeatureStatisticOutput t = new FeatureStatisticOutput(prefix+".maxPos");
		t.addCondition(new String[] {"transcript"}, "[U]");
		t.addCondition(new String[] {"Position"}, "[U]");
		program.add(t,"Position");
		
		new AlignedReadsDataToFeatureProgram(program).setProgress(context.getProgress(),()->"Identify maxpos").processStorage(reads);
		
		
		return null;
	}
	
	
}
