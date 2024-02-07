package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.function.Function;
import java.util.function.UnaryOperator;
import java.util.logging.Level;

import gedi.core.genomic.Genomic;
import gedi.grand3.Grand3Utils;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.PerrEstimator;
import gedi.grand3.estimation.estimators.FullModelEstimator;
import gedi.grand3.estimation.estimators.MaskErrorComponentModelEstimator;
import gedi.grand3.estimation.estimators.ModelEstimationMethod;
import gedi.grand3.estimation.estimators.ModelEstimator;
import gedi.grand3.estimation.models.Grand3Model;
import gedi.grand3.estimation.models.Grand3TruncatedBetaBinomialMixtureModel;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.KNMatrix;
import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.grand3.knmatrix.KnMatrixKey;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;
import gedi.util.sequence.Alphabet;

public class Grand3EstimateModel extends GediProgram {

	
	
	public Grand3EstimateModel(Grand3ParameterSet params) {
		addInput(params.nthreads);
		addInput(params.genomic);
		addInput(params.perrFile);
		addInput(params.knMatFile);
		addInput(params.subreadSemanticFile);
		addInput(params.subreadsToUse);
		addInput(params.experimentalDesignFile);
		addInput(params.noplot);
		addInput(params.confidence);
		addInput(params.profile);
		addInput(params.introntol);
		addInput(params.jointModel);
		addInput(params.estimMethod);
		
		addInput(params.prefix);
		addInput(params.test);
		addInput(params.debug);
		
		addOutput(params.modelFile);
		addOutput(params.modelBinFile);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		
		int pind = 0;
		int nthreads = getIntParameter(pind++);
		Genomic genomic = getParameter(pind++);
		File mmFile = getParameter(pind++); 
		File binomFile = getParameter(pind++); 
		File subreadFile = getParameter(pind++); 
		HashSet<String> useSubreads = new HashSet<>(getParameters(pind++));
		File designFile = getParameter(pind++); 
		boolean noplot = getBooleanParameter(pind++);
		boolean ci = getBooleanParameter(pind++);
		boolean profile = getBooleanParameter(pind++);
		int introntol = getIntParameter(pind++);
		String joinModels = getParameter(pind++);
		ArrayList<ModelEstimationMethod> estim = getParameters(pind++);

		
		String prefix = getParameter(pind++);
		String test = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		
		
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);
		String[] subreads = Grand3Utils.readSemantic(subreadFile);
		
		
		// which subreads to use
		boolean[] subreadsToUse = new boolean[subreads.length];
		if (useSubreads.size()==0)
			Arrays.fill(subreadsToUse, true);
		else
			for (int i=0; i<subreadsToUse.length; i++)
				subreadsToUse[i] = useSubreads.contains(i+"") || useSubreads.contains(subreads[i]);

		
		// subread,label,cond,min|max
		double[][][][] pre_perr = new double[subreads.length][design.getTypes().length][design.getNumSamples()][];

		PerrEstimator est = new PerrEstimator(mmFile, design, subreads);
		for (int t=0; t<design.getTypes().length; t++) {
			if (design.getSamplesNotHaving(design.getTypes()[t]).length>0)
				context.getLog().info(String.format("Estimate prior p.err for %s by regression using %s naive samples: %s",
						design.getTypes()[t].toString(),
						design.getTypes()[t].toString(),
						EI.wrap(design.getSamplesNotHaving(design.getTypes()[t])).map(i->design.getSampleNameForSampleIndex(i)).concat(",")));
			else
				context.getLog().info(String.format("Estimate prior p.err for %s without %s naive samples (median method).",
						design.getTypes()[t].toString(),design.getTypes()[t].toString()));
				
			for (int s=0; s<subreads.length; s++) 
				for (int i=0; i<design.getNumSamples(); i++)
					if (design.getLabelForSample(i,design.getTypes()[t])!=null)
						pre_perr[s][t][i] = est.estimate(i, s, design.getTypes()[t]);
		}
		
		
		HashMap<KnMatrixKey, KNMatrixElement[]> dataMap = Grand3Model.dataFromStatFile(binomFile.getPath());
		
		if (joinModels!=null) {
			if (joinModels.length()!=StringUtils.split(design.getSampleName(0),'.').length || !Alphabet.buildAlphabet("01").isValid(joinModels))
				throw new RuntimeException("Invalid specification of jointModel!");
			
			int[] keep = EI.seq(0,joinModels.length()).filterInt(i->joinModels.charAt(i)=='1').toIntArray();
			UnaryOperator<String> newName = s->{ String[] parts = StringUtils.split(s, '.'); return EI.wrap(keep).mapInt(i->parts[i]).concat(".");};
			context.getLog().info("Joining models (e.g.: "+design.getSampleName(0)+" -> "+newName.apply(design.getSampleName(0))+") ...");
			HashMap<KnMatrixKey, KNMatrixElement[]> dataMap2 = new HashMap<KnMatrixKey, KNMatrixElement[]>();
			for (KnMatrixKey k : dataMap.keySet()) {
				KnMatrixKey newKey = new KnMatrixKey(newName.apply(k.getCondition()), k.getSubread(), k.getLabel());
				dataMap2.put(newKey, dataMap2.containsKey(newKey) ? KNMatrix.sum(dataMap2.get(newKey),dataMap.get(k)) : dataMap.get(k));
			}
			for (KnMatrixKey k : new ArrayList<>(dataMap.keySet())) {
				KnMatrixKey newKey = new KnMatrixKey(newName.apply(k.getCondition()), k.getSubread(), k.getLabel());
				dataMap.put(k, dataMap2.get(newKey));
			}
		}
		
		Function<int[],KNMatrixElement[]> stiToData = sti->dataMap.get(new KnMatrixKey(design.getSampleNameForSampleIndex(sti[2]),sti[0],design.getTypes()[sti[1]].toString()));
		
		LineWriter out = getOutputWriter(0);
		boolean headerWritten = false;
		
		LineWriter profOut;
		if (profile) {
			profOut = new LineOrientedFile(prefix+".model.profile.tsv.gz").write();
			profOut.writef("Estimator\tCondition\tSubread\tLabel\tParameter\t%s\tdeltaLL\n",EI.wrap(new Grand3TruncatedBetaBinomialMixtureModel().createPar().getParamNames()).concat("\t"));
		}
		else profOut = null;
		
		context.getLog().info("Estimating models...");
		
		for (int mi=0; mi<estim.size(); mi++) {
			
			ModelEstimationMethod method = estim.get(mi);
			
			ModelEstimator estimator;
			if (method.equals(ModelEstimationMethod.Full))
				estimator=new FullModelEstimator(design,stiToData,pre_perr, subreads);
			else if (method.equals(ModelEstimationMethod.MaskErrorComponent))
				estimator=new MaskErrorComponentModelEstimator(design,stiToData,pre_perr, subreads);
			else
				throw new RuntimeException("Unknown estimation method: "+method);
			
			context.getLog().info("Estimating models using "+method+" method...");
			context.getLog().info(" with prior p.err distribution: "+estimator.getPriorPerr());
			
			ArrayList<ModelStructure> models = 
				EI.seq(0, subreads.length)
				.unfold(s->EI.seq(0,design.getTypes().length).map(t->new int[]{s,t}))
				.unfold(p->EI.seq(0, design.getNumSamples()).map(i->new int[]{p[0],p[1],i}))
				.filter(sti->design.getLabelForSample(sti[2], design.getTypes()[sti[1]])!=null)
				.parallelized(nthreads, 1, ei->ei.unfold(estimator::estimate))
				.sideEffect(
						m->!m.getBinom().isConverged(),
						m->context.warningf("Nelder-Mead did not converge on "+design.getTypes()[m.getType()].toString()+" binom model for "+design.getSampleNameForSampleIndex(m.getSample())+" ("+subreads[m.getSubread()]+")" )
						)
				.sideEffect(
						m->!m.getTBBinom().isConverged(),
						m->context.warningf("Nelder-Mead did not converge on "+design.getTypes()[m.getType()].toString()+" tbbinom model for "+design.getSampleNameForSampleIndex(m.getSample())+" ("+subreads[m.getSubread()]+")" )
						)
				.progress(context.getProgress(), design.getTypes().length*subreads.length*design.getNumSamples(), (r)->"Estimated "+design.getTypes()[r.getType()].toString()+" model for "+design.getSampleNameForSampleIndex(r.getSample())+" ("+subreads[r.getSubread()]+")")
				.removeNulls()
				.list();
			
			Collections.sort(models);
			for (ModelStructure m : models) {
				String status = m.getBinom().isConverged()?"OK":m.getBinom().getConvergence();
				if (!m.getBinom().isInBounds()) status = "OUT_OF_BOUNDS ("+status+")";
				context.logf("binom %s (%s %s): \t%s\t\t%s",design.getSampleNameForSampleIndex(m.getSample()),subreads[m.getSubread()],design.getTypes()[m.getType()].toString(),status,m.getBinom().toCreateString(true) );
			}
			for (ModelStructure m : models) {
				String status = m.getTBBinom().isConverged()?"OK":m.getTBBinom().getConvergence();
				if (!m.getTBBinom().isInBounds()) status = "OUT_OF_BOUNDS ("+status+")";
				context.logf("tbbinom %s (%s %s): \t%s\t\t%s",design.getSampleNameForSampleIndex(m.getSample()),subreads[m.getSubread()],design.getTypes()[m.getType()].toString(),status,m.getTBBinom().toCreateString(true) );
			}
			
			if (ci||profile) {
				context.getLog().info("Computing profile likelihood...");
				models = EI.wrap(models)
					.parallelized((int)(Math.ceil(nthreads)/8.0), 1, ei->ei.map(m->estimator.computeProfile(context.getLog(),profOut,m,nthreads,design,subreads)))
					.progress(context.getProgress(), models.size(), (r)->"Computed profile likelihood for the "+design.getTypes()[r.getType()].toString()+" model for "+design.getSampleNameForSampleIndex(r.getSample())+" ("+subreads[r.getSubread()]+")")
					.removeNulls()
					.list();
				Collections.sort(models);
			}
			
			if (!headerWritten) {
				models.get(0).writeTableHeader(out,ci||profile);
				headerWritten = true;
			}
			
			
			// if subread restriction, pay attention that the joint model is only estimated for these subreads
			// the paramteres for those are overwritten, and only those are printed into the table output
			// the other separate models are available in the binary output, but not used for target parameter estimation (see SubreadCounterKNMatrixPerTarget.count)
			if (estimator instanceof FullModelEstimator) {
			
				for (ModelStructure model : models) 
					model.writeTable(out,"Separate",ci||profile,design,subreads,pre_perr);
				
				int subreadsUsed = EI.seq(0,subreadsToUse.length).filterInt(i->subreadsToUse[i]).countInt();
				
				if (subreadsUsed>1) {
					context.getLog().info("Subreads used for joint model parameter estimation: "+EI.seq(0,subreads.length).filterInt(i->subreadsToUse[i]).map(i->subreads[i]).concat(","));
					context.getLog().info("Estimating joint models...");
					ArrayList<ModelStructure[]> blocks = EI.wrap(models)
							.filter(m->subreadsToUse[m.getSubread()])
							.multiplex(ModelStructure::compareNoS, ModelStructure.class).list();
	
					EI.wrap(blocks)
						.parallelized(nthreads, 1, ei->ei.map(block->{
							MaximumLikelihoodParametrization[] models2 = ((FullModelEstimator) estimator).estimateJoint(block, ci||profile);
							if (!models2[0].isConverged())
								context.warningf("Nelder-Mead did not converge on joint "+design.getTypes()[block[0].getType()].toString()+" binom model for "+design.getSampleNameForSampleIndex(block[0].getSample()) );
							if (!models2[1].isConverged())
								context.warningf("Nelder-Mead did not converge on joint "+design.getTypes()[block[0].getType()].toString()+" tbbinom model for "+design.getSampleNameForSampleIndex(block[0].getSample()) );
							return new int[] {block[0].getType(),block[0].getSample()};
						}))
						.progress(context.getProgress(), blocks.size(), (r)->"Estimated joint "+design.getTypes()[r[0]].toString()+" model for "+design.getSampleNameForSampleIndex(r[1]))
						.drain();
				}
				
				if (subreads.length>1) {
					for (ModelStructure model : models)
						// only models used for the joint model end up in the tabular output
						if (subreadsToUse[model.getSubread()])
							model.writeTable(out,"Joint",ci||profile,design,subreads,pre_perr);
				}
			}
			else if (estimator instanceof MaskErrorComponentModelEstimator) {
				for (ModelStructure model : models) 
					if (subreadsToUse[model.getSubread()])
						model.writeTable(out,"MaskedError",ci||profile,design,subreads,pre_perr);
				
			}
			
			if (mi==0) {
				PageFileWriter bin = new PageFileWriter(getOutputFile(1).getPath());
				for (ModelStructure model : models) 
					model.serialize(bin);
				bin.close();
			}
		}
		
		out.close();
		if (profile)
			profOut.close();
		
		
		
		if (!noplot) {
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".grand3.model.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/grand3.model.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		return null;
	}


	
}
