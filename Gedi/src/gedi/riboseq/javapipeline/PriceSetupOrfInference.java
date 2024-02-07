package gedi.riboseq.javapipeline;

import java.io.IOException;
import java.util.ArrayList;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionStorage;
import gedi.riboseq.inference.orf.OrfInference;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class PriceSetupOrfInference extends GediProgram {

	public PriceSetupOrfInference(PriceParameterSet params) {
		addInput(params.reads);
		addInput(params.genomic);
		addInput(params.introns);
		addInput(params.novelTranscripts);
		addInput(params.keepAnno);
		addInput(params.checkOrfs);
		addInput(params.checkAnnotation);
		
		addOutput(params.orfinference);
	}
	
	public String execute(GediProgramContext context) throws IOException {
		
		
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(0);
		Genomic genomic = getParameter(1);
		GenomicRegionStorage<Void> introns = getParameter(2);
		boolean novelTranscripts = getBooleanParameter(3);
		boolean keepAnno = getBooleanParameter(4);
		ArrayList<GenomicRegionStorage<Void>> checkOrfs = getParameters(5);
		boolean checkAnno = getBooleanParameter(6);
		
		OrfInference v = new OrfInference(genomic,reads);
		if (introns!=null && introns.size()>0)
			v.addSpliceJunctions(introns.ei());
		
		if (checkAnno) v.addCheckAnnotation();
		for (GenomicRegionStorage<Void> c : checkOrfs)
			v.addCheckOrfs(c.getName(), c.ei());
		
		v.setAllowNovelTranscripts(novelTranscripts);
		v.setRemoveAnno(!keepAnno);
		
		setOutput(0, v);
		
		return null;
	}
	
	

}
