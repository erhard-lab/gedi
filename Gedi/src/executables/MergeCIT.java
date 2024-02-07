package executables;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;

import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataMerger;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.reference.Chromosome;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.FileUtils;
import gedi.util.FunctorUtils;
import gedi.util.FunctorUtils.ParallellIterator;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.r.RRunner;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class MergeCIT {

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		if (args.length<3) {
			usage();
			System.exit(1);
		}
		
		boolean progress = false;
		boolean clear = false;
		HashSet<String> skip = new HashSet<>(); 
		String out = null;
		String referenceSequence = null;
		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else if (args[i].equals("-c"))
				clear = true;
			else if (args[i].equals("-ref"))
				referenceSequence = args[++i];
			else if (args[i].equals("-s"))
				EI.split(args[++i], ',').toCollection(skip);
			else {
				out = args[i++];
				args = Arrays.copyOfRange(args, i, args.length);
				i = args.length;
			}
		}
		if (out==null) {
			usage();
			System.exit(1);
		}
		
		EI.wrap(args).map(File::new).throwArg(File::exists,"File %s does not exist!");
		
		if (args.length==1) {
			if (new File(out).equals(new File(args[0])))
				System.exit(0);
			if (clear) {
				new File(args[0]).renameTo(new File(out));
				if (new File(args[0]+".metadata.json").exists())
					new File(args[0]+".metadata.json").renameTo(new File(out+".metadata.json"));
			}
			else {
				Files.copy(Paths.get(args[0]),Paths.get(out));
				if (new File(args[0]+".metadata.json").exists())
					Files.copy(Paths.get(args[0]+".metadata.json"),Paths.get(out+".metadata.json"));
			}
			
			System.exit(0);
		}
		
		CenteredDiskIntervalTreeStorage<AlignedReadsData>[] storages = (CenteredDiskIntervalTreeStorage[]) 
				EI.wrap(args)
					.map(CenteredDiskIntervalTreeStorage::load)
					.toArray(new CenteredDiskIntervalTreeStorage[0]);
		boolean compressed = EI.wrap(storages).mapToInt(cit->cit.isCompressed()?1:0).sum()>=storages.length/2;
		
		boolean barcodes = storages[0].getRandomRecord() instanceof BarcodedAlignedReadsData;
		Class<DefaultAlignedReadsData> cls = barcodes?(Class)BarcodedAlignedReadsData.class:(Class)DefaultAlignedReadsData.class;
		
		DynamicObject[] metas = EI.wrap(storages).map(CenteredDiskIntervalTreeStorage::getMetaData).filter(d->!d.isNull()).toArray(DynamicObject.class);
		DynamicObject meta = metas.length==storages.length?DynamicObject.merge(metas):DynamicObject.getEmpty();
		
		ExtendedIterator<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>>[] iterators;
		
		if (referenceSequence==null) 
			iterators = (ExtendedIterator[]) 
					EI.wrap(storages)
						.map(s->s.ei()
								.iff(skip.size()>0, ei->ei.filter(rgr->!skip.contains(rgr.getReference().getName())))
								.checkOrder((Comparator)FunctorUtils.naturalComparator(),()->s.getPath()))
						.toArray(new ExtendedIterator[0]);
		else {
			System.out.println("Will only merge on "+referenceSequence);
			String creferenceSequence = referenceSequence;
			iterators = (ExtendedIterator[]) 
			EI.wrap(storages)
				.map(s->s.ei(Chromosome.obtain(creferenceSequence,true),Chromosome.obtain(creferenceSequence,false))
						.iff(skip.size()>0, ei->ei.filter(rgr->!skip.contains(rgr.getReference().getName())))
						.checkOrder((Comparator)FunctorUtils.naturalComparator(),()->s.getPath()))
				.toArray(new ExtendedIterator[0]);
		}
		
		int[] conditions = EI.wrap(storages).mapToInt(s->s.getRandomRecord().getNumConditions()).toIntArray();
		AlignedReadsDataMerger merger = new AlignedReadsDataMerger(conditions);
		
		ParallellIterator<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>> pit = (ParallellIterator<ImmutableReferenceGenomicRegion<? extends AlignedReadsData>>) FunctorUtils.parallellIterator((Iterator[])iterators, FunctorUtils.naturalComparator(), ImmutableReferenceGenomicRegion.class);
		
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> outCit = new CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData>(out, cls,compressed);
		outCit.fill(pit.map(merger::merge).iff(progress, ei->ei.progress(new ConsoleProgress(System.err),-1,e->e.toLocationString()+e.getData())),new ConsoleProgress(System.err));
		if (!meta.isNull())
			outCit.setMetaData(meta);
		
		if (clear) {
			for (String a : args) {
				new File(a).delete();
				new File(a+".metadata.json").delete();
			}
		}
		
		long nbc = EI.wrap(storages).map(cit->FileUtils.getExtensionSibling(cit.getPath(), "barcodes.tsv")).filter(f->new File(f).exists()).count();
		if (nbc!=storages.length && nbc>0) Gedi.getLog().severe("Did find barcodes.tsv for part of the input files!");
		if (nbc==storages.length) {
			ExtendedIterator<String> it = EI.lines(FileUtils.getExtensionSibling(storages[0].getPath(), "barcodes.tsv"));
			for (int i=1; i<storages.length; i++)
				it = it.chain(EI.lines(FileUtils.getExtensionSibling(storages[i].getPath(), "barcodes.tsv")).skip(1));
			it.print(FileUtils.getExtensionSibling(outCit.getPath(), "barcodes.tsv"));
			
			if (clear) {
				EI.wrap(storages).map(cit->FileUtils.getExtensionSibling(cit.getPath(), "barcodes.tsv")).forEachRemaining(s->new File(s).delete());
			}
		}
		
		File[] readstats = EI.wrap(args).map(c->new File(new File(new File(c).getParentFile(),"report"),FileUtils.getNameWithoutExtension(c)+".reads.tsv")).toArray(File.class);
		if (EI.wrap(readstats).mapToInt(f->f.exists()?1:0).sum()==args.length) {
			String name = FileUtils.getNameWithoutExtension(out);
			String script = new File(new File(new File(out).getParentFile(),"report"),name+".reads.R").getPath();
			out = new File(new File(new File(out).getParentFile(),"report"), name+".reads.tsv").getPath();
			
			RRunner r = new RRunner(script);
			r.setNumeric("files", "c("+EI.wrap(readstats).map(f->"\""+f.getPath()+"\"").concat(",")+")");
			r.set("out", out);
			r.addSource(r.getClass().getResourceAsStream("/resources/R/mergeread.R"));
			r.run(false);
			
			String tsv = out;
			for (String png : EI.fileNames(new File(new File(out).getParentFile().getParentFile(),"report").getPath()).filter(f->f.startsWith(name+".reads") && f.endsWith(".png")).loop()) {
				String title = FileUtils.getNameWithoutExtension(png);
				FileUtils.writeAllText("{\"plots\":[{\"section\":\"Mapping summary\",\"id\":\"ID"+png.replace('.', '_')+"\",\"title\":\""+title+"\",\"description\":\"Mapping summary. For a description see mapping statistics!\",\"img\":\""+png+"\",\"script\":\""+script+"\",\"csv\":\""+tsv+"\"}]}",new File("report/"+FileUtils.getExtensionSibling(png, ".report.json")));
			}

		}
		
		
	}

	private static void usage() {
		System.out.println("MergeCIT [-c] [-p] [-s skip1,skip2,...] <output> <file1> <file2> ... \n\n -c removes the input files after successful merging\n -p shows progress\n -s skip chromosomes");
	}
	
}
