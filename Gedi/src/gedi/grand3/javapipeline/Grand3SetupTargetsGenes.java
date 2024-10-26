package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.grand3.targets.TargetCollection;
import gedi.grand3.targets.geneExonic.GeneExonicTargetCollection;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3SetupTargetsGenes extends GediProgram {

	public static enum ExonIntronMode {
		Standard{
			public boolean getExonForGlobal() {return true;	}
			public boolean getIntronForGlobal() {return false;}
			public boolean getExonForTarget() {return true;}
			public boolean getIntronForTarget() {return false;}
		},
		Talkseq {
			public boolean getExonForGlobal() {return false;	}
			public boolean getIntronForGlobal() {return true;}
			public boolean getExonForTarget() {return true;}
			public boolean getIntronForTarget() {return true;}
		}, 
		Intronic {
			public boolean getExonForGlobal() {return false;	}
			public boolean getIntronForGlobal() {return true;}
			public boolean getExonForTarget() {return false;}
			public boolean getIntronForTarget() {return true;}
		};
		
		public abstract boolean getExonForGlobal();
		public abstract boolean getIntronForGlobal();
		public abstract boolean getExonForTarget();
		public abstract boolean getIntronForTarget();
		
	}
	
	
	public Grand3SetupTargetsGenes(Grand3ParameterSet params) {
		
		addInput(params.genomic);
		addInput(params.introntol);
		addInput(params.mode);
		addInput(params.overlap);
		addInput(params.genemode);
//		addInput(params.targetMergeTable);
		
		setRunFlag(params.targets, "genes");
		
		addOutput(params.targetCollection);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int pind = 0;
		Genomic genomic = getParameter(pind++);
		int introntol = getIntParameter(pind++);
		ReadCountMode mode = getParameter(pind++);
		ReadCountMode overlap = getParameter(pind++);
		ExonIntronMode exin = getParameter(pind++);
		String mergetab = getParameter(pind++);
		
		context.getLog().info("Parameters will be estimated for genes ("+exin+")!");
		
		
		boolean exonForGlobal = exin.getExonForGlobal();
		boolean intronForGlobal = exin.getIntronForGlobal();
		boolean exonForOutput = exin.getExonForTarget();
		boolean intronForOutput = exin.getIntronForTarget();
		
		String first = genomic.getOriginList().get(0);
		TargetCollection targets =  new GeneExonicTargetCollection(
				genomic,
				genomic.getGenes(),
				r->r.isMitochondrial()?"Mito":genomic.getOrigin(r).getId(),
				r->genomic.getLength(r),
				s->s.equals(first),
				exonForGlobal,intronForGlobal,exonForOutput,intronForOutput,
				genomic.getTranscripts(),
				mode,
				overlap,
				introntol);
		
		targets.checkValid();
		
//		if (mergetab!=null) {
//			HeaderLine h = new HeaderLine();
//			HashMap<String, String> map = EI.lines(mergetab).header(h).split('\t').index(a->a[h.apply("name")], a->a[h.apply("merged")]);
//			context.getLog().info("Loaded merge table with "+map.size()+" entries.");
//			targets = new TargetCollectionMappedName(targets, a->map.getOrDefault(a, a));
//		}
		
		setOutput(0, targets);

		
		return null;
	}


	
}
