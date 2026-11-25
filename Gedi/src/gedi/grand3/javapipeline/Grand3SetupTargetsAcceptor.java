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

public class Grand3SetupTargetsAcceptor extends GediProgram {

	


	public Grand3SetupTargetsAcceptor(Grand3ParameterSet params) {
		
		addInput(params.genomic);
		addInput(params.introntol);
		addInput(params.mode);
		addInput(params.overlap);
		addInput(params.minOverlap);
		addInput(params.genemode);
		
		setRunFlag(params.targets, "acceptor");
		
		addOutput(params.targetCollection);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		Genomic genomic = getParameter(pind++);
		int introntol = getIntParameter(pind++);
		ReadCountMode mode = getParameter(pind++);
		ReadCountMode overlap = getParameter(pind++);
		int minOverlap = getIntParameter(pind++);
		ExonIntronMode exin = getParameter(pind++);
		
		context.getLog().info("Parameters will be estimated for splice acceptor sites!");
		
		String first = genomic.getOriginList().get(0);
		TargetCollection targets =  new SpliceAcceptorTargetCollection(
				genomic,
				genomic.getGenes(),
				r->r.isMitochondrial()?"Mito":genomic.getOrigin(r).getId(),
				r->genomic.getLength(r),
				s->s.equals(first),
				exin.getExonForGlobal(),exin.getIntronForGlobal(),exin.getExonForTarget(),exin.getIntronForTarget(),
				genomic.getTranscripts(),
				mode,
				overlap,
				introntol,
				minOverlap);
		
		targets.checkValid();
		
		setOutput(0, targets);
		
		return null;
	}
	
}
