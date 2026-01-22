package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.grand3.targets.TargetCollection;
import gedi.grand3.targets.geneExonic.AllExonIntronTargetCollection;
import gedi.grand3.targets.geneExonic.GeneExonicTargetCollection;
import gedi.grand3.targets.geneExonic.SpliceAcceptorTargetCollection;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3SetupTargetsExonIntron extends GediProgram {

	
	

	public Grand3SetupTargetsExonIntron(Grand3ParameterSet params) {
		
		addInput(params.genomic);
		addInput(params.introntol);
		addInput(params.mode);
		addInput(params.overlap);
		
		setRunFlag(params.targets, "exin");
		
		addOutput(params.targetCollection);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		Genomic genomic = getParameter(pind++);
		int introntol = getIntParameter(pind++);
		ReadCountMode mode = getParameter(pind++);
		ReadCountMode overlap = getParameter(pind++);
		
		context.getLog().info("Parameters will be estimated for exons and introns!");
		
		String first = genomic.getOriginList().get(0);
		TargetCollection targets =  new AllExonIntronTargetCollection(
				genomic,
				genomic.getGenes(),
				r->r.isMitochondrial()?"Mito":genomic.getOrigin(r).getId(),
				r->genomic.getLength(r),
				s->s.equals(first),
				false,true,true,true,
				genomic.getTranscripts(),
				mode,
				overlap,
				introntol);
		
		targets.checkValid();
		
		setOutput(0, targets);
		
		return null;
	}


	
}
