package gedi.grand3;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.estimators.ModelEstimationMethod;
import gedi.grand3.estimation.models.Grand3BinomialMixtureModel;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.targets.TargetCollection;
import gedi.grand3.targets.geneExonic.GeneExonicTargetCollection;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public class Grand3Utils {

		
	public static Strandness getStrandness(File file) throws IOException {
		return Strandness.valueOf(EI.lines(file).split('\t').filter(a->a[0].equals("Used")).getUniqueResult(true, true)[1]);
	}

	
	public static TargetCollection getTargets(Genomic genomic, Strandness strandness, ReadCountMode mode, ReadCountMode overlap, int introntol) {
		String first = genomic.getOriginList().get(0);
		return new GeneExonicTargetCollection(
				genomic,
				genomic.getGenes(),
				r->r.isMitochondrial()?"Mito":genomic.getOrigin(r).getId(),
				r->genomic.getLength(r),
				s->s.equals(first),
				true,false,true,false,
				genomic.getTranscripts(),
				mode,
				overlap,
				introntol);
	}



	public static void writeSemantic(String[] subreads, File f) throws IOException {
		LineWriter wr = new LineOrientedFile(f.getPath()).write();
		wr.writeLine("Subread\tSemantic");
		for (int s=0; s<subreads.length; s++)
			wr.writef("%d\t%s\n", s, subreads[s]);
		wr.close();
	}
	
	public static String[] readSemantic(File f) throws IOException {
		ArrayList<String> re = new ArrayList<String>();
		for (String[] a : EI.lines(f).split('\t').skip(1).loop()) {
			if (Integer.parseInt(a[0])!=re.size())throw new RuntimeException("SubreadSemantic is not sequential!");
			re.add(a[1]);
		}
		return re.toArray(new String[0]);
	}
	
	public static ModelStructure[][][] readModelsTsv(File modelFile, ExperimentalDesign design, String[] subreads, ModelEstimationMethod method) throws IOException {
		
		String methname = method.equals(ModelEstimationMethod.Full)?"Joint":"MaskedError";
		
		
		HashMap<String,Integer> sampleIndex = new HashMap<String, Integer>();
		HashMap<String,Integer> subreadIndex = EI.wrap(subreads).indexPosition();
		HashMap<String,Integer> typeIndex = new HashMap<String, Integer>();
		
		for (int i=0; i<design.getNumSamples(); i++)
			sampleIndex.put(design.getSampleNameForSampleIndex(i), i);
		for (int i=0; i<design.getTypes().length; i++)
			typeIndex.put(design.getTypes()[i].toString(), i);
		
		Grand3BinomialMixtureModel binom = new Grand3BinomialMixtureModel();
		Grand3TruncatedBetaBinomialMixtureModel tbbinom = new Grand3TruncatedBetaBinomialMixtureModel();
		
		int maxs = 0;
		int maxt = 0;
		ArrayList<ModelStructure> l = new ArrayList<ModelStructure>();
		
		HeaderLine h = new HeaderLine();
		for (String[] a : EI.lines(modelFile).header(h).split('\t').loop()) {
			if (a[h.get("Estimator")].equals(methname)) {
				ModelStructure m = new ModelStructure(binom.createPar(), tbbinom.createPar());
				
				m.parse(a,h,sampleIndex,subreadIndex,typeIndex);
				maxs = Math.max(maxs, m.getSubread());
				maxt = Math.max(maxt, m.getType());
				l.add(m);
			}
		}
		
		
		ModelStructure[][][] re = new ModelStructure[maxs+1][maxt+1][design.getNumSamples()];
		for (ModelStructure m : l)
			re[m.getSubread()][m.getType()][m.getSample()] = m;
		
		return re;
	}
	

	public static ModelStructure[][][] readModels(File modelFile, int numCond) throws IOException {
		PageFile f = new PageFile(modelFile.getPath());
		
		Grand3BinomialMixtureModel binom = new Grand3BinomialMixtureModel();
		Grand3TruncatedBetaBinomialMixtureModel tbbinom = new Grand3TruncatedBetaBinomialMixtureModel();
		
		int maxs = 0;
		int maxt = 0;
		ArrayList<ModelStructure> l = new ArrayList<ModelStructure>();
		while (!f.eof()) {
			ModelStructure m = new ModelStructure(binom.createPar(), tbbinom.createPar());
			m.deserialize(f);
			maxs = Math.max(maxs, m.getSubread());
			maxt = Math.max(maxt, m.getType());
			l.add(m);
		}
		
		ModelStructure[][][] re = new ModelStructure[maxs+1][maxt+1][numCond];
		for (ModelStructure m : l)
			re[m.getSubread()][m.getType()][m.getSample()] = m;
		
		return re;
	}
	
	
	
}
