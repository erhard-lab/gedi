package gedi.grand3.javapipeline;


import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.function.Function;
import java.util.function.UnaryOperator;

import gedi.app.Gedi;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.grand3.Grand3Utils;
import gedi.grand3.estimation.TargetEstimationResult;
import gedi.grand3.estimation.TargetEstimationResult.ModelType;
import gedi.grand3.estimation.TargetEstimationResult.SingleEstimationResult;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.grand3.experiment.PseudobulkDefinition;
import gedi.grand3.targets.TargetCollection;
import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableInteger;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

public class Grand3OutputFlatfiles extends GediProgram {

	
	
	public Grand3OutputFlatfiles(Grand3ParameterSet params, boolean hasMappedTarget) {
		addInput(params.genomic);
		if (hasMappedTarget)
			addInput(params.pseudobulkBinFile);
		else
			addInput(params.targetBinFile);
		addInput(params.experimentalDesignFile);
		addInput(params.prefix);
		addInput(params.targetCollection);
		addInput(params.debug);
		addInput(params.forceDense);
		addInput(params.forceSparse);
		addInput(params.strandnessFile);
		addInput(params.pseudobulkFile);
		addInput(params.pseudobulkName);
		addInput(params.pseudobulkMinimalPurity);
		addInput(params.targetMergeTable);
		addInput(params.outputMixBeta);

		
		if (hasMappedTarget)
			addOutput(params.pseudobulkFolder);
		else
			addOutput(params.targetFolder);
	}
	
	
	
	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		
		int pind = 0;
		Genomic genomic = getParameter(pind++);
		File targetFile = getParameter(pind++); 
		File designFile = getParameter(pind++); 
		
		String prefix = getParameter(pind++);
		TargetCollection targets = getParameter(pind++);
		boolean debug = getBooleanParameter(pind++);
		boolean forceDense = getBooleanParameter(pind++);
		boolean forceSparse = getBooleanParameter(pind++);
		File strandnessFile = getParameter(pind++); 
		String pseudobulkFile = getParameter(pind++);
		String pseudobulkName = getParameter(pind++);
		double pseudobulkMinimalPurity = getDoubleParameter(pind++);
		String targetMergeTab = getParameter(pind++);
		boolean writeMix = getBooleanParameter(pind++);
		
		
		Strandness strandness = Grand3Utils.getStrandness(strandnessFile);
		ExperimentalDesign design = ExperimentalDesign.fromTable(designFile);
		
		
		File folder = getOutputFile(0);
		folder.mkdirs();
		if (!folder.isDirectory()) throw new IOException("Could not create directory "+folder);
		
		
		String[] columnNames;
		int[] columnToSample;
		if (pseudobulkFile!=null) {
			PseudobulkDefinition psdef = new PseudobulkDefinition(pseudobulkFile,design,null,pseudobulkMinimalPurity);
			columnNames = psdef.getPseudobulkNames();
			columnToSample=psdef.getSampleForPseudobulks();
		}
		else {
			columnNames = EI.seq(0, design.getCount()).map(i->design.getFullName(i)).toArray(String.class);
			columnToSample = design.getIndexToSampleId();
		}
		
		boolean sparse = columnNames.length>30;
		if (forceDense) sparse=false;
		if (forceSparse) sparse=true;
		
		HashMap<String, ArrayList<String>> map = new HashMap<String, ArrayList<String>>();
		if (targetMergeTab!=null) {
			HeaderLine h = new HeaderLine();
			map = EI.lines(targetMergeTab).header(h).split('\t').indexMulti(a->a[h.apply("merged")], a->a[h.apply("name")]);
		}
		
		if (sparse) 
			outputSparse(context,genomic,design,columnNames,columnToSample,targetFile,folder,writeMix);
		else
			outputDense(context,genomic,design,columnNames,columnToSample,targetFile,folder,map,writeMix);

		return null;
	} 
	
	private void outputSparse(GediProgramContext context, Genomic genomic, ExperimentalDesign design, String[] columnNames, int[] columnToSample, File targetFile, File dir, boolean writeMix) throws IOException {
		
		HashMap<String,MutableInteger> g2l = new HashMap<>();
		genomic.getTranscripts().ei().forEachRemaining(tr->g2l.computeIfAbsent(tr.getData().getGeneId(), x->new MutableInteger()).max(tr.getRegion().getTotalLength())); 

		
		Function<String, String> sym = genomic.getGeneTable("geneId", "symbol");
		context.getLog().info("Writing sparse matrices...");

		PageFile pf = new PageFile(targetFile.getPath());
		int storedNumFeatures = pf.getInt();
		long storedNumCounts = pf.getLong();
		long[] storedNumNtr = new long[design.getTypes().length];
		for (int i=0; i<storedNumNtr.length; i++)
			storedNumNtr[i] = pf.getLong();
		
		EI.seq(0,columnNames.length).map(i->columnNames[i]).print(new File(dir,"barcodes.tsv.gz").getPath());
		LineWriter features = new LineOrientedFile(dir, "features.tsv.gz").write();
		LineWriter counts = new LineOrientedFile(dir,"matrix.mtx.gz").write();
		counts.writeLine("%%MatrixMarket matrix coordinate integer general");
		counts.writeLine("%metadata_json: {\"format_version\": 2, \"software_version\": \"grand3_"+Gedi.version("GRAND3")+"\"}");
		counts.writef("%d %d %d\n",storedNumFeatures,columnNames.length,storedNumCounts);
		LineWriter[][] ntrs = mkWriter(dir,"ntr","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][] ntrLower = mkWriter(dir,"lower","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][] ntrUpper = mkWriter(dir,"upper","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][] ntrAlpha = mkWriter(dir,"alpha","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][] ntrBeta = mkWriter(dir,"beta","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][] ntrInte = mkWriter(dir,"integral","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[] shape = mkWriter2(dir,"shape","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[] shapeLLR = mkWriter2(dir,"llr","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[] shapeLL = mkWriter2(dir,"ll","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr);
		LineWriter[][][] mix = writeMix?new LineWriter[][][]{
			mkWriter(dir,"mix","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr),
			mkWriter(dir,"alpha1","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr),
			mkWriter(dir,"beta1","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr),
			mkWriter(dir,"alpha2","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr),
			mkWriter(dir,"beta2","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr),
			mkWriter(dir,"bmintegral","real",design.getTypes(),storedNumFeatures,columnNames.length,storedNumNtr)
		}:null;
		
		int numFeatures = 0;
		long numCounts = 0;
		long[] numNtr = new long[design.getTypes().length];
		
		for (TargetEstimationResult res : pf.ei(()->new TargetEstimationResult())
				.progress(context.getProgress(),storedNumFeatures,t->t.getTarget().toString())
				.loop()) {
		if (!res.isExonic()) System.out.println("Sparse output for introns not implemented yet!");
		
			numFeatures++; // must start with 1 in mm file, so this is correct
			int len = 0;
			if (g2l.get(res.getTarget())!=null)
				len = g2l.get(res.getTarget()).N;
			features.writef("%s\t%s\tGene Expression\t%s\t%d\n",res.getTarget(),sym.apply(res.getTarget()),res.getGenome(),len);
			
			IntIterator it = res.iterateCounts();
			while (it.hasNext()) {
				int i = it.nextInt();
				int barcode = i+1;
				counts.writef("%d %d %.0f\n", numFeatures,barcode,res.getCountOrZero(i));
				numCounts++;
			}
			
			for (int t=0; t<design.getTypes().length; t++) {
				it = res.iterate(t);
				while (it.hasNext()) {
					int i = it.nextInt();
					int barcode = i+1;
					SingleEstimationResult val = res.get(t,i);
					
					for(ModelType type : ModelType.values()) {
						ntrs[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getMap());
						ntrLower[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getLower());
						ntrUpper[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getUpper());
						ntrAlpha[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getAlpha());
						ntrBeta[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBeta());
						ntrInte[t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getIntegral());
						if (writeMix) {
							mix[0][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixMix());
							mix[1][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixAlpha1());
							mix[2][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixBeta1());
							mix[3][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixAlpha2());
							mix[4][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixBeta2());
							mix[5][t][type.ordinal()].writef("%d %d %.4f\n", numFeatures,barcode,val.getModel(type).getBetaMixIntegral());
						}
					}
					shape[t].writef("%d %d %.4f\n", numFeatures,barcode,val.getShape().getShape());
					shapeLLR[t].writef("%d %d %.4f\n", numFeatures,barcode,val.getShape().getLogLikShape()-val.getShape().getLogLikGlobal());
					shapeLL[t].writef("%d %d %.4f\n", numFeatures,barcode,val.getShape().getLogLikShape());
					numNtr[t]++;
				}
			}
		}
		pf.close();
		
		features.close();
		finishWriter(counts);
		finishWriters(ntrs);
		finishWriters(ntrLower);
		finishWriters(ntrUpper);
		finishWriters(ntrAlpha);
		finishWriters(ntrBeta);
		finishWriters(ntrInte);
		finishWriters(shape);
		finishWriters(shapeLLR);
		finishWriters(shapeLL);
		if (writeMix) 
			for (int x=0; x<mix.length; x++)
				finishWriters(mix[x]);
		
		if (storedNumCounts!=numCounts) throw new RuntimeException("Targets file corrupted: NumCounts does not match!");
		for (int i=0; i<storedNumNtr.length; i++)
			if (storedNumNtr[i]!=numNtr[i]) throw new RuntimeException("Targets file corrupted: NumNTR does not match!");
	}

	private void finishWriters(LineWriter[] wr) throws IOException {
		for (int t=0; t<wr.length; t++) finishWriter(wr[t]);
	}
	private void finishWriters(LineWriter[][] wr) throws IOException {
		for (int t=0; t<wr.length; t++) for (LineWriter w : wr[t]) finishWriter(w);
	}
	private void finishWriter(LineWriter wr) throws IOException {
		wr.close();
	}


	private LineWriter[][] mkWriter(File folder, String name, String format, MetabolicLabelType[] types, int nx, int ny, long[] nz) throws IOException {
		LineWriter[][] re = new LineWriter[types.length][ModelType.values().length];
		for (int t=0; t<types.length; t++)
			for (int m=0; m<ModelType.values().length; m++) {
				re[t][m] = new LineOrientedFile(folder, types[t].toString()+"."+ModelType.values()[m].name()+"."+name+".mtx.gz").write();
				re[t][m].writeLine("%%MatrixMarket matrix coordinate "+format+" general");
				re[t][m].writeLine("%metadata_json: {\"format_version\": 2, \"software_version\": \"grand3_"+Gedi.version("GRAND3")+"\"}");
				re[t][m].writef("%d %d %d\n",nx,ny,nz[t]);
			}
		return re;
	}
	
	private LineWriter[] mkWriter2(File folder, String name, String format, MetabolicLabelType[] types, int nx, int ny, long[] nz) throws IOException {
		LineWriter[] re = new LineWriter[types.length];
		for (int t=0; t<types.length; t++) {
				re[t] = new LineOrientedFile(folder, types[t].toString()+"."+name+".mtx.gz").write();
				re[t].writeLine("%%MatrixMarket matrix coordinate "+format+" general");
				re[t].writeLine("%metadata_json: {\"format_version\": 2, \"software_version\": \"grand3_"+Gedi.version("GRAND3")+"\"}");
				re[t].writef("%d %d %d\n",nx,ny,nz[t]);
			}
		return re;
	}



	private void outputDense(GediProgramContext context, Genomic g, ExperimentalDesign design, String[] columnNames, int[] columnToSample, File targetFile, File folder, HashMap<String, ArrayList<String>> map, boolean writeMix) throws IOException {
		context.getLog().info("Writing flat files (dense)...");
		
		HashMap<String,MutableInteger> e2l = new HashMap<>();
		g.getTranscripts().ei().forEachRemaining(tr->e2l.computeIfAbsent(tr.getData().getGeneId(), x->new MutableInteger()).max(tr.getRegion().getTotalLength())); 
		for (String merged : map.keySet())
			e2l.put(merged, new MutableInteger(EI.wrap(map.get(merged)).mapToInt(ge->e2l.getOrDefault(ge, new MutableInteger()).N).sum()));
			
		HashMap<String,MutableInteger> i2l = new HashMap<>();
		g.getTranscripts().ei().forEachRemaining(tr->i2l.computeIfAbsent(tr.getData().getGeneId(), x->new MutableInteger()).max(tr.getRegion().invert().getTotalLength())); 
		for (String merged : map.keySet())
			i2l.put(merged, new MutableInteger(EI.wrap(map.get(merged)).mapToInt(ge->i2l.getOrDefault(ge, new MutableInteger()).N).sum()));
		
		
		PageFile pf = new PageFile(targetFile.getPath());
		
		int storedNumFeatures = pf.getInt();
		long storedNumCounts = pf.getLong();
		long[] storedNumNtr = new long[design.getTypes().length];
		for (int i=0; i<storedNumNtr.length; i++)
			storedNumNtr[i] = pf.getLong();
		
		LineOrientedFile exonicFile = new LineOrientedFile(folder, "data.tsv.gz");
		LineOrientedFile intronicFile = new LineOrientedFile(folder, "intronic.tsv.gz");
		
		LineWriter exonic = exonicFile.write(); 
		LineWriter intronic = intronicFile.write(); 
		writeDenseHeader(exonic,design,columnNames,columnToSample,writeMix);
		writeDenseHeader(intronic,design,columnNames,columnToSample,writeMix);
		boolean hadExonic = false;
		boolean hadIntronic = false;
		
		for (TargetEstimationResult res : pf.ei(()->new TargetEstimationResult()).loop()) {
			if (res.isEmpty()) continue;
			
			if (res.isExonic()) hadExonic = true;
			else if (res.isIntronic()) hadIntronic = true;
			else throw new RuntimeException("Invalid category!");
			
			LineWriter out = res.isExonic()?exonic:intronic;
			HashMap<String,MutableInteger> lenmap = res.isExonic()?e2l:i2l;
			
			String sym = g.getGeneTable("symbol").apply(res.getTarget());
			if (sym==null) {
				if (res.getTarget().contains("_")) {
					sym = g.getGeneTable("symbol").apply(res.getTarget().substring(0,res.getTarget().indexOf(("_"))));
					if (sym!=null)
						sym = sym+res.getTarget().substring(res.getTarget().indexOf(("_")));	
				}
				
				if (sym==null)
					sym = res.getTarget();
					
			}
			out.writef("%s\t%s\t%s\t%d",
					res.getTarget(),
					sym,
					res.getGenome(),lenmap.getOrDefault(res.getTarget(),new MutableInteger()).N
					);
			for (int i=0; i<columnNames.length; i++)
				out.writef("\t%.1f",res.getCountOrZero(i));
			for (ModelType type : ModelType.values())
				for (int t=0; t<design.getTypes().length; t++)
					for (int i=0; i<columnNames.length; i++)
						if (design.getLabelForSample(columnToSample[i], design.getTypes()[t])!=null) {
							if (res.get(t,i)!=null) {
								out.writef("\t%.3f\t%.3f\t%.3f\t%.4f\t%.4f\t%.4f",
										res.get(t,i).getModel(type).getLower(),
										res.get(t,i).getModel(type).getMap(),
										res.get(t,i).getModel(type).getUpper(),
										res.get(t,i).getModel(type).getAlpha(),
										res.get(t,i).getModel(type).getBeta(),
										res.get(t,i).getModel(type).getIntegral()
										);
							}
							else
								out.writef("\tNA\tNA\tNA\tNA\tNA\tNA");
							if (writeMix) {
								if (res.get(t,i)!=null) {
									out.writef("\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
											res.get(t,i).getModel(type).getBetaMixMix(),
											res.get(t,i).getModel(type).getBetaMixAlpha1(),
											res.get(t,i).getModel(type).getBetaMixBeta1(),
											res.get(t,i).getModel(type).getBetaMixAlpha2(),
											res.get(t,i).getModel(type).getBetaMixBeta2(),
											res.get(t,i).getModel(type).getBetaMixIntegral()
											);
								}
								else
									out.writef("\tNA\tNA\tNA\tNA\tNA");
							}
						}
			for (int t=0; t<design.getTypes().length; t++)
				for (int i=0; i<columnNames.length; i++)
					if (design.getLabelForSample(columnToSample[i], design.getTypes()[t])!=null) {
						if (res.get(t,i)!=null) 
							out.writef("\t%.2f\t%.1f\t%.1f",
									res.get(t,i).getShape().getShape(), 
									res.get(t,i).getShape().getLogLikShape()-res.get(t,i).getShape().getLogLikGlobal(),
									res.get(t,i).getShape().getLogLikShape());
						else
							out.writef("\tNA\tNA\tNA");
					} 
			out.writeLine();
			
		}
		pf.close();
		
		exonic.close();
		intronic.close();
		
		if (!hadExonic) exonicFile.delete();
		if (!hadIntronic) intronicFile.delete();
		
		if (!hadExonic && hadIntronic) intronicFile.renameTo(exonicFile);
	}



	private void writeDenseHeader(LineWriter out, ExperimentalDesign design, String[] columnNames, int[] columnToSample, boolean writeMix) throws IOException {
		out.writef("Gene\tSymbol\tCategory\tLength");
		for (int i=0; i<columnNames.length; i++)
			out.writef("\t%s Read count",columnNames[i]);
		for (ModelType type : ModelType.values())
			for (MetabolicLabelType t: design.getTypes())
				for (int i=0; i<columnNames.length; i++)
					if (design.getLabelForSample(columnToSample[i], t)!=null) {
						out.writef("\t%s %s %s NTR Lower CI\t%s %s %s NTR MAP\t%s %s %s NTR Upper CI\t%s %s %s alpha\t%s %s %s beta\t%s %s %s integral",
								columnNames[i],t.toString(),type.toString(),
								columnNames[i],t.toString(),type.toString(),
								columnNames[i],t.toString(),type.toString(),
								columnNames[i],t.toString(),type.toString(),
								columnNames[i],t.toString(),type.toString(),
								columnNames[i],t.toString(),type.toString());
						if (writeMix)
							out.writef("%s %s %s mix\t%s %s %s alpha1\t%s %s %s beta1\t%s %s %s alpha2\t%s %s %s beta2\t%s %s %s bmintegral",
									columnNames[i],t.toString(),type.toString(),
									columnNames[i],t.toString(),type.toString(),
									columnNames[i],t.toString(),type.toString(),
									columnNames[i],t.toString(),type.toString(),
									columnNames[i],t.toString(),type.toString(),
									columnNames[i],t.toString(),type.toString());
					}
		for (MetabolicLabelType t: design.getTypes())
			for (int i=0; i<columnNames.length; i++)
				if (design.getLabelForSample(columnToSample[i], t)!=null)
					out.writef("\t%s Shape\t%s LLR\t%s LL",columnNames[i],columnNames[i],columnNames[i]);
		out.writeLine();
	}



	private void writeDenseGene(LineWriter out, Genomic genomic, ExperimentalDesign design, TargetEstimationResult[] results) throws IOException {
		// either all nor nothing!
		int nulls = 0;
		for (TargetEstimationResult r : results) if (r!=null) nulls++;
		
		if (nulls==0) return;
		if (nulls!=results.length) throw new RuntimeException("Data file corrupted: Not all types are available!");
		
		
		
	}


	
}
