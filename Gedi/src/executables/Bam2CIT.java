package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import gedi.app.Gedi;
import gedi.app.extension.ExtensionContext;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStorageCapabilities;
import gedi.core.region.GenomicRegionStorageExtensionPoint;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.region.bam.BamGenomicRegionStorage;
import gedi.region.bam.BamGenomicRegionStorage.PairedEndHandling;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.ParseUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableInteger;
import gedi.util.userInteraction.progress.ConsoleProgress;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


public class Bam2CIT {

	
	private static int checkMultiParam(String[] args, int index, ArrayList<String> re)  {
		while (index<args.length && !args[index].startsWith("-")) 
			re.add(args[index++]);
		return index-1;
	}

	private static int checkIntParam(String[] args, int index) {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new RuntimeException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(args[index]);
	}
	private static String checkParam(String[] args, int index)  {
		if (index>=args.length || args[index].startsWith("-")) throw new RuntimeException("Missing argument for "+args[index-1]);
		return args[index];
	}
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		if (args.length<2) {
			usage();
			System.exit(1);
		}
		
		Gedi.startup(false);
		
		
		boolean is10x = false;
		boolean isDropseq = false;
		boolean isUmi = false;
		boolean isUmiONT = false;
		String umiPattern = null;
		boolean umiAllowMulti = false;
		boolean progress = false;
		boolean keepIds = false;
		boolean compress = false;
		boolean var = true;
		boolean keepMito = false;
		boolean nosec = false;
		boolean unspec = false;
		boolean anti = false;
		boolean join = false;
		String name = null;
		String removePref = null;
		HashMap<String,String> barcodeList = null;
		int minmaq = -1;
		int sechip = -1;
		int head = -1;
		Genomic check = null;
		String out = null;
		PairedEndHandling peh = null;

		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else if (args[i].equals("-check")) {
				ArrayList<String> gnames = new ArrayList<>();
				i = checkMultiParam(args, ++i, gnames);
				if (gnames.size()==0) throw new RuntimeException("No genomic given!");
				check = Genomic.get(gnames);
			} else if (args[i].equals("-id"))
				keepIds = true;
			else if (args[i].equals("-head"))
				head = checkIntParam(args,++i);
			else if (args[i].equals("-pe"))
				peh = ParseUtils.parseEnumNameByPrefix(checkParam(args,++i),true,PairedEndHandling.class);
			else if (args[i].equals("-compress"))
				compress = true;
			else if (args[i].equals("-novar"))
				var = false;
			else if (args[i].equals("-name"))
				name = checkParam(args, ++i);
			else if (args[i].equals("-minmaq"))
				minmaq = checkIntParam(args, ++i);
			else if (args[i].equals("-sechip"))
				sechip = checkIntParam(args, ++i);
			else if (args[i].equals("-nosec"))
				nosec = true;
			else if (args[i].equals("-strandunspecific"))
				unspec = true;
			else if (args[i].equals("-anti"))
					anti = true;
			else if (args[i].equals("-")){}
			else if (args[i].equals("-barcodelist")) {
				HeaderLine h = new HeaderLine();
				barcodeList = EI.lines(checkParam(args, ++i)).header(h).split("\t").index(a->StringUtils.removeFooter(a[h.get("Barcode")],"-1"),a->a[h.get("Sample")]);
			} 
			else if (args[i].equals("-removePrefix")) {
				removePref = checkParam(args, ++i);
			} 
			else if (args[i].equals("-10x")) {
				is10x = true;
			} 
			else if (args[i].equals("-dropseq")) {
				isDropseq = true;
			} 
			else if (args[i].equals("-umi")) 
				isUmi = true;
			else if (args[i].equals("-umiOnt")) 
				isUmiONT = true;
			else if (args[i].equals("-umiAllowMulti")) {
				umiAllowMulti = true;
			} else if (args[i].equals("-umiPattern")) {
				isUmi = true;
				umiPattern = checkParam(args, ++i);
			} else if (args[i].equals("-keepMito"))
				keepMito = true;
			else if (args[i].equals("-join"))
				join = true;
			else if (args[i].startsWith("-"))
				throw new IllegalArgumentException("Parameter "+args[i]+" unknown!");
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
		if (!out.endsWith(".cit")) out = out+".cit";
		
		EI.wrap(args).map(File::new).throwArg(File::exists,"File %s does not exist!");
		args = EI.wrap(args).unfold(s->s.endsWith(".bamlist")?EI.lines2(s,"#"):EI.singleton(s).str()).toArray(String.class); 
		
		BamGenomicRegionStorage storage = new BamGenomicRegionStorage(args);
		int numCond = storage.getRandomRecord().getNumConditions();

		
		if (anti)
			storage.setStrandness(Strandness.Antisense);
		if (peh!=null)
			storage.setPairedEndHandling(peh);
		
		LineWriter incons = null;  
		if (check!=null)
			storage.check(check,incons = new LineOrientedFile(out+".inconsistent").write());
		storage.setIgnoreVariations(!var);
		storage.setOnlyPrimary(nosec);
		if (join) 
			storage.setJoinMates(true);
		
		storage.setKeepReadNames(keepIds);
		Class<?> dataClass = DefaultAlignedReadsData.class;

		if (isUmiONT) {
			storage.setUmiOnt();
			dataClass = BarcodedAlignedReadsData.class;
		}
		
		if (isUmi) {
			if (umiPattern!=null)
				storage.setUmi(umiPattern);
			else
				storage.setUmi();
			if (umiAllowMulti)
				storage.setMinimalAlignmentQuality(0);
			dataClass = BarcodedAlignedReadsData.class;
		}
		
		if (is10x) {
			dataClass = BarcodedAlignedReadsData.class;
			
			File[] bcs = EI.wrap(args)
				.map(bam->{
					File f= FileUtils.getExtensionSibling(new File(bam), "filtered_barcodes.tsv");
					if (!f.exists()) 
						f= new File(new File(new File(bam).getAbsoluteFile().getParentFile(),"filtered_feature_bc_matrix"), "barcodes.tsv.gz");
					if (!f.exists()) 
						f= new File(new File(new File(bam).getAbsoluteFile().getParentFile(),"sample_filtered_feature_bc_matrix"), "barcodes.tsv.gz");
					if (!f.exists()) 
						f= new File(new File(new File(new File(new File(bam).getAbsoluteFile().getParentFile(),"Solo.out"),"Gene"),"filtered"), "barcodes.tsv");
					if (f.exists())
						Gedi.getLog().info("Found filtered barcodes for "+bam+" in "+f);
					else
						Gedi.getLog().warning("Did not find filtered barcodes for "+bam+"!");
					return f;
				}).toArray(File.class);
			
			if (new File(FileUtils.getExtensionSibling(out, ".barcodes.tsv")).exists()) {
				Gedi.getLog().warning("Will not create barcodes file, file already present!");
			}
			else if (EI.wrap(bcs).filter(f->f.exists()).count()==bcs.length) {
				Gedi.getLog().info("Creating barcodes file!");
				HashMap<String,String> nBarcodeList = new HashMap<String, String>();
				String[] conds = storage.getMetaDataConditions();
				if (conds.length>1) throw new RuntimeException("Do not call with multiple bams, use MergeCIT instead!");
				if (name!=null) conds[0] = name;
				try (LineWriter wr = new LineOrientedFile(FileUtils.getExtensionSibling(out, ".barcodes.tsv")).write()) {
					wr.writeLine("Library\tBarcode\tSample");
					for (int c=0; c<conds.length; c++) {
						for (String bc : EI.lines(bcs[c]).map(bc->StringUtils.removeFooter(bc, "-1")).loop())
							if (barcodeList!=null) {
								if (barcodeList.containsKey(bc))
									wr.writef("%s\t%s\t%s\n",conds[c],bc,barcodeList.get(bc));
							} else {
								wr.writef("%s\t%s\t%s\n",conds[c],bc,conds[c]);
								nBarcodeList.put(bc, conds[c]);
							}
						
					}
				}
				if (barcodeList==null) barcodeList = nBarcodeList;
				
			}
			else {
				Gedi.getLog().warning("Will not create barcodes file, not all filtered_feature_bc_matrix found!");
			}
			storage.set10x(barcodeList);
			keepMito = true;
		}
		if (isDropseq) {
			dataClass = BarcodedAlignedReadsData.class;
			
			if (minmaq<0) minmaq = 255;
			if (args.length!=1) throw new RuntimeException("Only a single bam supported right now!");
			if (barcodeList==null) barcodeList = determineBarcodes(args[0],minmaq);
			storage.setDropseq(barcodeList);
			
			System.err.println("Creating barcodes file!");
		
			String[] conds = storage.getMetaDataConditions();
			try (LineWriter wr = new LineOrientedFile(FileUtils.getExtensionSibling(out, ".barcodes.tsv")).write()) {
				wr.writeLine("Library\tBarcode\tSample");
				for (String bc : EI.wrap(barcodeList.keySet()).loop())
					wr.writef("%s\t%s\t%s\n",conds[0],bc,barcodeList.get(bc));
			}
			keepMito = true;
		}

		if (minmaq>=0)
			storage.setMinimalAlignmentQuality(minmaq);
		
		@SuppressWarnings("rawtypes")
		GenomicRegionStorage outStorage = GenomicRegionStorageExtensionPoint.getInstance().get(new ExtensionContext().add(Boolean.class, compress).add(String.class, out).add(Class.class, dataClass), GenomicRegionStorageCapabilities.Disk, GenomicRegionStorageCapabilities.Fill);
		NumericArray mitocount = NumericArray.createMemory(numCond, NumericArrayType.Double);
		
		if (head>0 || !keepMito || sechip>0 || unspec || removePref!=null) {
			ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> it = null;
			
			String uRemovePref = removePref;
			if (removePref==null) it = storage.ei();
			else it = EI.wrap(storage.getReferenceSequences()).
					filter(r->r.getName().startsWith(uRemovePref)).unfold(r->storage.ei(r)).
					map(r->new ImmutableReferenceGenomicRegion<>(gedi.core.reference.Chromosome.obtain(r.getReference().getName().substring(uRemovePref.length()),r.getReference().getStrand()), r.getRegion(), r.getData()));

			
			if (head>0) it = it.head(head);
			if (!keepMito) it = it.filter(r->{
				boolean mito = r.getReference().isMitochondrial();
				if (mito) r.getData().addTotalCountsForConditions(mitocount, ReadCountMode.Weight);
				return !mito;
			});
			
			if (sechip>0) {
				int usechip = sechip;
				it=it.map(read->new ImmutableReferenceGenomicRegion<>(read.getReference(), read.mapMaybeOutSide(new ArrayGenomicRegion(0,usechip)),read.getData()));
			}
			
			if (unspec)
				it=it.map(read->read.toStrandIndependent());
			
			if (progress) it = it.progress(new ConsoleProgress(System.err),-1,r->r.toLocationString());
			outStorage.fill(it);
		} else {
			outStorage.fill(storage,progress?new ConsoleProgress(System.err):null);
		}
		
		if (incons!=null)
			incons.close();
		
		DynamicObject meta = storage.getMetaData();
//		meta = meta.cascade(DynamicObject.from("conditions", 
//				DynamicObject.from(new Object[] {
//						DynamicObject.from("readlengths", storage.getReadLengths())
//				})
//				));
		if (storage.getNumConditions()==1 && name!=null)
			meta = meta.cascade(DynamicObject.from("conditions", 
					DynamicObject.from(new Object[] {
							DynamicObject.from("name",name)
					})
					));
		else if (storage.getNumConditions()==1 && meta.get(".conditions[0].name").isNull())
			meta = meta.cascade(DynamicObject.from("conditions", 
					DynamicObject.from(new Object[] {
							DynamicObject.from("name",storage.getName())
					})
					));
		outStorage.setMetaData(meta);
		
		if (!keepMito) {
			System.out.println("Condition\tMitochondrial");
			String[] conds = storage.getMetaDataConditions();
			for (int c=0; c<conds.length; c++) {
				System.out.print(conds[c]);
				System.out.println("\t"+mitocount.getDouble(c));
			}
			System.out.println();
		}
		System.out.println("--------------------");
		System.out.println("Total\t"+mitocount.sum());
		
	}

	private static HashMap<String, String> determineBarcodes(String sam, int minmaq) throws IOException {
		SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(sam));
		HashMap<String, MutableInteger> counter = new HashMap<String, MutableInteger>();
		for (SAMRecord r : EI.wrap(reader.iterator()).progress(r->"Counting cell barcodes: "+r.getReferenceName()+":"+r.getAlignmentStart()+" Barcodes="+counter.size()).loop())
			if (!r.getReadUnmappedFlag() && r.getMappingQuality()>=minmaq) {
				String bc = r.getStringAttribute("XC");
				counter.computeIfAbsent(bc, x->new MutableInteger()).N++;
			}
		reader.close();
		
		double[] counts = new double[counter.size()];
		String[] bcs = new String[counter.size()];
		int index = 0;
		for (String bc : counter.keySet()) {
			counts[index]=counter.get(bc).N;
			bcs[index++]=bc;
		}
		ArrayUtils.parallelSort(counts, bcs);
		ArrayUtils.reverse(counts);
		ArrayUtils.reverse(bcs);
		int ind = kneedle(counts, -1);
		HashMap<String,String> re = new HashMap<String, String>();
		long n=0;
		for (int i=0; i<ind; i++) {
			re.put(bcs[i], "");
			n+=counter.get(bcs[i]).N;
		}
		System.err.printf("Found %d barcodes with a total of %d UMIs, median %.0f UMIs!\n",ind,n,counts[ind/2]);
		return re;
	}

	private static double dist2d(double[] a, double[] b, double[] c) {
		double[] v1 = ArrayUtils.subtract(b.clone(), c);
		double[] v2 = ArrayUtils.subtract(a.clone(), b);
		double det = v1[0]*v2[1]-v1[1]*v2[0];
		double s  = v1[0]*v1[0]+v1[1]*v1[1];
		return det/Math.sqrt(s);
	}

	public static int kneedle(double[] values, int sign) {
		double[] start = {1,values[0]};
		double[] end = {values.length, values[values.length-1]};
		int re = 0;
		double reval = Double.NEGATIVE_INFINITY;
		for (int idx=1; idx<=values.length; idx++) {
			double test = sign*-1*dist2d(new double[] {idx,values[idx-1]},start,end);
			if (test>reval) {
				re = idx-1;
				reval = test;
			}
		}
		return re;
	}
	private static void usage() {
		System.out.println("Bam2CIT [-p] [-id] [-compress] [-minmaq <MAQ>] [-keepMito] [-novar] [-nosec] [-10x] [-umi [-umiAllowMulti] [-umiPattern <regex-all-groups-are-used>]] [-barcodelist <multiseq-table>] [-removePrefix <prefix>] <output> <file1> <file2> ... \n\n -p shows progress\n -id add ids to CIT\n -removePrefix filters reads for that and removes the prefix (e.g. for 10x runs with human/mouse combined)\n -barcodelist <multiseq-table>  needs to be a tsv file with columns Barcode and Sample!");
	}
	
}
