package gedi.macoco.javapipeline;

import java.util.HashMap;
import java.util.function.BiConsumer;

import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.math.stat.inference.isoforms.EquivalenceClassCountEM;
import gedi.util.program.GediProgramContext;

public class EquivalenceClassesMaximumLikelihood extends EstimateTpm {



	public EquivalenceClassesMaximumLikelihood(MacocoParameterSet params) {
		super(params,params.emMLTable);
	}



	@Override
	protected void estimate(GediProgramContext context, Genomic genomic, MemoryIntervalTreeStorage<Transcript> trans, Strandness strand, int cond, String condName, double[] eff, String[][] E, double[] counts,
			BiConsumer<String, Double> transUnnorm) {
		
		HashMap<String, Double> traToEffl = trans.ei().index(r->r.getData().getTranscriptId(), r->eff[r.getRegion().getTotalLength()]);

		EquivalenceClassCountEM<String> algo = new EquivalenceClassCountEM<String>(E, counts, traToEffl::get);
		algo.compute(2000, 10000, (t,a)->{
			transUnnorm.accept(t, a/traToEffl.get(t));
		});
		
	}
}
