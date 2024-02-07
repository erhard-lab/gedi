package gedi.grand3.javapipeline;


import java.io.IOException;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.grand3.targets.TargetCollection;
import gedi.grand3.targets.geneExonic.GeneExonicTargetCollection;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3SetupTargetsTest extends GediProgram {

	
	
	public Grand3SetupTargetsTest(Grand3ParameterSet params) {
		
		addInput(params.genomic);
		addInput(params.introntol);
		addInput(params.mode);
		addInput(params.overlap);
		addInput(params.test);
		
		setRunFlag(params.targets, "test");
		
		addOutput(params.targetCollection);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		Genomic genomic = getParameter(pind++);
		int introntol = getIntParameter(pind++);
		ReadCountMode mode = getParameter(pind++);
		ReadCountMode overlap = getParameter(pind++);
		String test = getParameter(pind++);

		
		
		context.getLog().info("Parameters will be estimated for "+test+"!");
		
		MemoryIntervalTreeStorage<String> genes = new MemoryIntervalTreeStorage<String>(String.class);
		MemoryIntervalTreeStorage<Transcript> tr = new MemoryIntervalTreeStorage<Transcript>(Transcript.class);
		ImmutableReferenceGenomicRegion<String> gg = ImmutableReferenceGenomicRegion.parse(genomic, test, "Test");
		genes.fill(genomic.getGenes().ei(gg));
		tr.fill(genomic.getTranscripts().ei(gg));
		
		TargetCollection targets =  new GeneExonicTargetCollection(
				genes,
				r->"Test",
				r->genomic.getLength(r),
				s->true,
				true,false,true,false,
				tr,
				mode,
				overlap,
				introntol) ;
		
		
		setOutput(0, targets);

		
		return null;
	}


	
}
