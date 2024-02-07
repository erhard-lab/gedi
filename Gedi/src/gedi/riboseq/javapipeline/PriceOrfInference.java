package gedi.riboseq.javapipeline;

import java.io.IOException;
import java.util.logging.Level;

import gedi.riboseq.inference.orf.NoiseModel;
import gedi.riboseq.inference.orf.OrfInference;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.riboseq.inference.orf.StartCodonScorePredictor;
import gedi.riboseq.utils.RiboUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.collections.doublecollections.DoubleIterator;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class PriceOrfInference extends GediProgram {

	public PriceOrfInference(PriceParameterSet params) {
		addInput(params.prefix);
		addInput(params.nthreads);
		addInput(params.orfinference);
		addInput(params.codons);
		addInput(params.startmodel);
		addInput(params.noisemodel);
		
		addOutput(params.pvals);
		addOutput(params.orfstsv);
		addOutput(params.orfsbin);
		
	}
	
	public String execute(GediProgramContext context) throws IOException {
		
		String prefix = getParameter(0);
		int nthreads = getIntParameter(1);
		OrfInference v = getParameter(2);
		int chunk = 10;
		
		
		
		
		PageFile fm = new PageFile(prefix+".start.model");
		StartCodonScorePredictor predictor = new StartCodonScorePredictor(); 
		predictor.deserialize(fm);
		fm.close();
		v.setStartCodonPredictor(predictor);
		
		
		fm = new PageFile(prefix+".noise.model");
		NoiseModel m = new NoiseModel();
		m.deserialize(fm);
		
		v.setNoiseModel(m);
		fm.close();
		
		
		
		String[] conditions = v.getConditions();

		
		context.getLog().log(Level.INFO, "Infer ORFs");
		LineWriter tab = new LineOrientedFile(prefix+".orfs.tsv").write();
		PriceOrf.writeTableHeader(tab, conditions);
		PageFileWriter tmp = new PageFileWriter(prefix+".orfs.bin");
		DoubleArrayList pvals = new DoubleArrayList();
		RiboUtils.processCodonsSink(prefix+".codons.bin", "ORF inference", context.getProgress(), ()->"Cache: "+StringUtils.toString(v.getNoiseModel().getCacheSize()),nthreads, chunk, PriceOrf.class, 
				ei->ei
					//.filter(r->r.getReference().getName().equals("1") || 
					//		r.intersects(ImmutableReferenceGenomicRegion.parse("16-:18807284-18807423|18808731-18808751")))
					.demultiplex(o->v.inferOrfs(o).ei()),
				n->{
					try {
//							if (n.getData().getExpP()>v.getTestThreshold() && n.getData().getAbortiveP()>v.getTestThreshold()) {
								n.toMutable().serialize(tmp);
								pvals.add(n.getData().getCombinedP());
//							}
							n.getData().writeTableLine(tab, n);
					} catch (IOException e) {
						throw new RuntimeException("Could not write ORFs!",e);
					}
					
				});
		tab.close();
		tmp.close();
		
		PageFileWriter pvf = new PageFileWriter(getOutputFile(0).getAbsolutePath()); 
		DoubleIterator dit = pvals.iterator();
		while( dit.hasNext()) pvf.putDouble(dit.nextDouble());
		pvf.close();
		
		
		return null;
	}
	

}
