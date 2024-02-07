package gedi.riboseq.javapipeline;

import java.io.IOException;
import java.util.function.Predicate;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.riboseq.inference.clustering.RiboClusterBuilder;
import gedi.riboseq.utils.RiboUtils;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class PriceClusterReads extends GediProgram {

	public PriceClusterReads(PriceParameterSet params) {
		addInput(params.prefix);
		addInput(params.reads);
		addInput(params.filter);
		addInput(params.genomic);
		addInput(params.nthreads);
		addOutput(params.clusters);
	}
	
	public String execute(GediProgramContext context) throws IOException {
		
		String prefix = getParameter(0);
		GenomicRegionStorage<AlignedReadsData> reads = getParameter(1);
		String filter = getParameter(2);
		Genomic genomic = getParameter(3);
		int nthreads = getIntParameter(4);
		
		Predicate<ReferenceGenomicRegion<AlignedReadsData>> ffilter = RiboUtils.parseReadFilter(filter);
		
		RiboClusterBuilder clb = new RiboClusterBuilder(prefix, reads, ffilter, genomic.getTranscripts(), 1, 5, context.getProgress(), nthreads);
		clb.setContinueMode(true);
		clb.setPredicate(r->genomic.getSequenceNames().contains(r.getReference().getName()));
		clb.build();
		
		return null;
	}
	
	

}
