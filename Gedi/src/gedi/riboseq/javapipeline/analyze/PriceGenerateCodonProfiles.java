package gedi.riboseq.javapipeline.analyze;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.zip.GZIPOutputStream;

import gedi.app.extension.ExtensionContext;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.riboseq.analysis.MajorIsoform;
import gedi.riboseq.analysis.PriceAnalysis;
import gedi.riboseq.analysis.PriceAnalysisExtensionPoint;
import gedi.riboseq.codonprocessor.RDataOutput;
import gedi.riboseq.inference.orf.PriceOrf;
import gedi.riboseq.javapipeline.PriceParameterSet;
import gedi.util.FileUtils;
import gedi.util.ReflectionUtils;
import gedi.util.datastructure.array.SparseMemoryFloatArray;
import gedi.util.datastructure.dataframe.DataFrame;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.ParallelizedIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.StringLineWriter;
import gedi.util.mutable.MutablePair;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.r.RDataWriter;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class PriceGenerateCodonProfiles extends GediProgram {

	
	public PriceGenerateCodonProfiles(PriceParameterSet params) {
		addInput(params.majorIsoformCit);
		addInput(params.nthreads);
		addInput(params.indices);
		addInput(params.genomic);
		addInput(params.orfs);
		addInput(params.prefix);
		
		setRunFlag(params.genProfile);
		
		addOutput(params.profileFile);
	}
	

	public String execute(GediProgramContext context) throws IOException, InterruptedException {
		int p=0;
		File mifile = getParameter(p++);
		int nthreads = getIntParameter(p++);
		File codonsFile = getParameter(p++);
		Genomic g = getParameter(p++);
		File orffile = getParameter(p++);
		String prefix = getParameter(p++);
		
		CenteredDiskIntervalTreeStorage<MajorIsoform> cit = new CenteredDiskIntervalTreeStorage<>(mifile.getAbsolutePath());
		CenteredDiskIntervalTreeStorage<PriceOrf> orfs = new CenteredDiskIntervalTreeStorage<>(orffile.getAbsolutePath());
		CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codons = new CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray>(codonsFile.getAbsolutePath());
		
		cit.ei().progress(context.getProgress(), (int)cit.size(), r->r.toLocationString())
			.writeRDS(getOutputFile(0).getPath(), (int)cit.size(), nthreads,5,mi->mi.getData().getTranscript().getData().getGeneId(), (mi,out)->write(mi,out,codons,orfs,g));
		
		return null;
	}

	private static void write(ImmutableReferenceGenomicRegion<MajorIsoform> mi, RDataWriter writer, CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codons, CenteredDiskIntervalTreeStorage<PriceOrf> orfs, Genomic genomic) throws IOException {
		ImmutableReferenceGenomicRegion<Transcript> trans = mi.getData().getTranscript();
		double[][] cmat = new double[trans.getRegion().getTotalLength()][];
		int ncond = -1;
		for (ImmutableReferenceGenomicRegion<SparseMemoryFloatArray> codon : codons.ei(trans).loop()) {
			if (trans.getRegion().containsUnspliced(codon.getRegion())) {
				ArrayGenomicRegion cod = trans.induce(codon.getRegion());
				cmat[cod.getStart()] = codon.getData().toDoubleArray();
				ncond = cmat[cod.getStart()].length;
			}
		}
		if (ncond==-1)
			ncond = codons.getRandomRecord().length();
		for (int i=0; i<cmat.length; i++) if(cmat[i]==null) cmat[i] = new double[ncond];
		
		LinkedList<MutablePair<GenomicRegion,PriceOrf>> horfs = orfs.ei(trans)
					.filter(o->trans.getRegion().containsUnspliced(o.getRegion()))
					.map(o->new MutablePair<GenomicRegion,PriceOrf>(trans.induce(o.getRegion()),o.getData()))
					.toList();
		
		
		String seq = genomic.getSequence(trans).toString();
		
		String gene = trans.getData().getGeneId();
		String transcript = trans.getData().getTranscriptId();
		String symbol = genomic.getGeneTable("geneId", "symbol").apply(trans.getData().getGeneId());
		
		writer.startList(null, 4);
		writer.write(null, new String[] {gene,symbol,transcript,trans.toLocationString()});
		writer.write(null, new String[] {seq});
		writer.write(null, cmat);
		{
			writer.startList(null, horfs.size());
			for (MutablePair<GenomicRegion,PriceOrf> o : horfs) {
				
				writer.startList(null, 3);
				writer.write(null,new int[] {o.Item1.getStart()});
				writer.write(null,new int[] {o.Item1.getEnd()});
				writer.write(null,new double[] {o.Item2.getCombinedP()});
				writer.endList("start","end","pval");
			}
			writer.endList(EI.wrap(horfs).map(o->o.Item2.getTableId()).toArray(String.class));	
		}
		
		writer.endList("ids","sequence","codons","orfs");
		
	}
	
	
	public static void main(String[] args) throws IOException {
		
//		
//		RDataWriter out = new RDataWriter(new GZIPOutputStream(new FileOutputStream("test.rds")));
//		out.writeHeader(false);
//		out.startList(null, 6);
//		out.write(null, new int[] {5,3});
//		
//		out.startList(null, 2);
//		out.write(null, new String[] {"x"});
//		out.write(null, new String[] {"y"});
//		out.endList(new String[] {"X","Y"});
//		
//		out.write(null, new String[] {"a"});
//		out.write(null, new boolean[] {true});
//		out.writeNull();
//		out.write(null, new double[] {2.1});
//		out.endList(new String[] {"a","b","c","d","e","f"});
//		out.finish();
		
//		RDataWriter out = new RDataWriter(new GZIPOutputStream(new FileOutputStream("test2.rds")));
//		out.writeHeader(false);
//		out.startList(null, 2);
//		out.write(null, new double[] {0});
//		out.write(null, new double[] {1});
//		out.endList(new String[] {"E0","E1"});
//		out.finish();
//		
//		
//		int n = 2;
//		
//		out = new RDataWriter(new GZIPOutputStream(new FileOutputStream("test.rds")));
//		out.writeHeader(false);
//		out.startList(null, n);
//		
//		ExtendedIterator<MutablePair<String,byte[]>> pit = EI.seq(0, 2).map(i->{
//			try {
//				ByteArrayOutputStream stream = new ByteArrayOutputStream(1<<16);
//				RDataWriter hout = new RDataWriter(stream);
//				hout.write(null, new double[] {i});
//				stream.flush();
//				return new MutablePair<>("E"+i,stream.toByteArray());
//			} catch (IOException e) {
//				throw new RuntimeException(e);
//			}
//			
//
//		});
//		ArrayList<String> names = new ArrayList<String>();
//		for (MutablePair<String,byte[]> mi : pit.loop()) {
//			names.add(mi.Item1);
//			out.writeRaw(mi.Item2);
//		}
//		
//		out.endList(names.toArray(new String[0]));
//		out.finish();
//		
//		EI.seq(0, 200).writeRDS("test3.rds",200,10,1,i->"E"+i,(i,o)->{o.write(null, new double[] {i});Thread.sleep(new Random().nextInt(300));});
		
		
		CenteredDiskIntervalTreeStorage<PriceOrf> orfs = new CenteredDiskIntervalTreeStorage<>("price/stressors.orfs.cit");
		CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray> codons = new CenteredDiskIntervalTreeStorage<SparseMemoryFloatArray>("price/stressors.codons.cit");
		CenteredDiskIntervalTreeStorage<MajorIsoform> cit = new CenteredDiskIntervalTreeStorage<>("price/stressors.majorisoform.cit");
		Genomic g = Genomic.get("h.ens90");
		int n = 150;
		cit.ei().head(n).progress(new ConsoleProgress(), (int)n, r->r.toLocationString())
		.writeRDS("test.rds", (int)n,5,5, mi->mi.getData().getTranscript().getData().getGeneId(), (mi,out)->{
			write(mi, out, codons, orfs, g);
		});
	
	}

}

