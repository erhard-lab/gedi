package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.region.bam.BamGenomicRegionStorage;
import gedi.util.FileUtils;
import gedi.util.StringUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IterateIntoSink;
import gedi.util.io.text.HeaderLine;
import gedi.util.sequence.DnaSequence;
import gedi.util.userInteraction.progress.ConsoleProgress;

/**
 * Bam2CIT for Parse Biosciences (Evercode, STARsolo CB_UMI_Complex) data that is split into
 * several sublibraries (one BAM per sublibrary).
 * <p>
 * All BAMs are converted into a CIT with a single condition. Since the same cell barcode
 * (round1_round2_round3) occurs in several sublibraries, each cell barcode gets a two nucleotide
 * prefix that identifies the sublibrary:
 * sublibrary 1 -> AA, 2 -> AC, 3 -> AG, 4 -> AT, 5 -> CA, ... 16 -> TT
 * (i.e. sublibrary k is k-1 written in base 4 with digits ACGT).
 * <p>
 * The barcode stored in the CIT is: prefix (2nt) + CB without underscores (24nt) + UB (10nt).
 * <p>
 * Only reads with a CB contained in the called cells table (auto_cell_calling, columns
 * sublibrary, CB_hex, CB_dT) of the very same sublibrary are kept. The sublibrary of a BAM file
 * is inferred from its path (default pattern: sublibrary_(\d+)).
 * <p>
 * With -split, several CITs are written in a single pass over the BAMs: one per called cells
 * table X.called_cells.tsv (written to X.cit). A cell may be contained in several tables.
 * <p>
 * The matching barcodes.tsv (Library/Barcode/Sample, next to the CIT) must contain the prefixed
 * barcodes, i.e. make_atlas_tables.R has to prepend the same prefix.
 */
public class ParseSeqBam2CIT {

	private static final char[] NT = {'A','C','G','T'};

	public static String sublibraryPrefix(int sublibrary) {
		if (sublibrary<1 || sublibrary>16) throw new IllegalArgumentException("Sublibrary must be in 1..16: "+sublibrary);
		int k = sublibrary-1;
		return new String(new char[] {NT[k/4],NT[k%4]});
	}

	private static String checkParam(String[] args, int index)  {
		if (index>=args.length || args[index].startsWith("-")) throw new RuntimeException("Missing argument for "+args[index-1]);
		return args[index];
	}
	private static int checkIntParam(String[] args, int index) {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new RuntimeException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(re);
	}

	/**
	 * One output CIT with its called cells.
	 */
	private static class Group {
		String cellsFile;
		String out;
		String name;
		HashSet<String>[] whitelist; // per bam file: CBs without underscores
		HashSet<String> prefixedCalled = new HashSet<>();
		AtomicLong kept = new AtomicLong();
		AtomicLong droppedNotCalled = new AtomicLong();
		AtomicLong droppedLength = new AtomicLong();
		CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData> storage;
		IterateIntoSink<ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData>> sink;
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		if (args.length<2) {
			usage();
			System.exit(1);
		}

		Gedi.startup(false);

		boolean progress = false;
		boolean keepIds = false;
		boolean var = true;
		boolean nosec = false;
		boolean anti = false;
		int minmaq = -1;
		int head = -1;
		String name = null;
		String cellsFile = null;
		String[] splitFiles = null;
		String[] cbColumns = {"CB_hex","CB_dT"};
		Pattern sublibPattern = Pattern.compile("sublibrary_(\\d+)");
		String out = null;

		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else if (args[i].equals("-id"))
				keepIds = true;
			else if (args[i].equals("-novar"))
				var = false;
			else if (args[i].equals("-nosec"))
				nosec = true;
			else if (args[i].equals("-anti"))
				anti = true;
			else if (args[i].equals("-minmaq"))
				minmaq = checkIntParam(args, ++i);
			else if (args[i].equals("-head"))
				head = checkIntParam(args, ++i);
			else if (args[i].equals("-name"))
				name = checkParam(args, ++i);
			else if (args[i].equals("-cells"))
				cellsFile = checkParam(args, ++i);
			else if (args[i].equals("-split"))
				splitFiles = StringUtils.split(checkParam(args, ++i), ',');
			else if (args[i].equals("-cbColumns"))
				cbColumns = StringUtils.split(checkParam(args, ++i), ',');
			else if (args[i].equals("-sublibraryPattern"))
				sublibPattern = Pattern.compile(checkParam(args, ++i));
			else if (args[i].equals("-")){}
			else if (args[i].startsWith("-"))
				throw new IllegalArgumentException("Parameter "+args[i]+" unknown!");
			else {
				if (splitFiles==null) out = args[i++];
				args = Arrays.copyOfRange(args, i, args.length);
				i = args.length;
			}
		}

		// output groups
		ArrayList<Group> groups = new ArrayList<>();
		if (splitFiles!=null) {
			if (cellsFile!=null) throw new IllegalArgumentException("Do not use -cells together with -split!");
			for (String f : splitFiles) {
				if (f.endsWith(".bam") || f.endsWith(".bamlist"))
					throw new IllegalArgumentException("-split got a bam file ("+f+") instead of called cells tables; is the list of called_cells.tsv files empty?");
				Group g = new Group();
				g.cellsFile = f;
				g.out = f.endsWith(".called_cells.tsv") ? f.substring(0, f.length()-".called_cells.tsv".length())+".cit" : FileUtils.getExtensionSibling(f, ".cit");
				groups.add(g);
			}
		}
		else {
			if (out==null || cellsFile==null) {
				usage();
				System.exit(1);
			}
			if (!out.endsWith(".cit")) out = out+".cit";
			Group g = new Group();
			g.cellsFile = cellsFile;
			g.out = out;
			groups.add(g);
		}
		if (args.length==0) {
			usage();
			System.exit(1);
		}
		HashSet<String> outs = new HashSet<>();
		for (Group g : groups) {
			if (!outs.add(g.out)) throw new IllegalArgumentException("Output "+g.out+" given more than once!");
			for (String a : args)
				if (new File(a).getAbsoluteFile().equals(new File(g.out).getAbsoluteFile()))
					throw new IllegalArgumentException("Output "+g.out+" is also an input!");
		}

		EI.wrap(args).map(File::new).throwArg(File::exists,"File %s does not exist!");
		args = EI.wrap(args).unfold(s->s.endsWith(".bamlist")?EI.lines2(s,"#"):EI.singleton(s).str()).toArray(String.class);

		// sublibrary per bam file
		int[] sublib = new int[args.length];
		HashSet<Integer> seenSublib = new HashSet<>();
		for (int f=0; f<args.length; f++) {
			Matcher m = sublibPattern.matcher(new File(args[f]).getAbsolutePath());
			int s = -1;
			while (m.find()) s = Integer.parseInt(m.group(1));
			if (s<0) throw new RuntimeException("Cannot infer sublibrary from path "+args[f]+" (pattern "+sublibPattern.pattern()+")");
			if (!seenSublib.add(s)) throw new RuntimeException("Sublibrary "+s+" given more than once!");
			sublib[f] = s;
			Gedi.getLog().info("Sublibrary "+s+" (prefix "+sublibraryPrefix(s)+"): "+args[f]);
		}
		String[] prefixes = new String[args.length];
		for (int f=0; f<args.length; f++) prefixes[f] = sublibraryPrefix(sublib[f]);

		// called cells per group; union with underscores (as in the bam) for the prefilter in the storage
		HashMap<String,String> unionWithUnderscores = new HashMap<>();
		int cbLength = -1;
		for (Group g : groups) {
			HashMap<Integer,HashSet<String>> cellsPerSublib = new HashMap<>();
			HeaderLine h = new HeaderLine();
			ExtendedIterator<String[]> cit = EI.lines(g.cellsFile).header(h).split("\t");
			ArrayList<String> needed = new ArrayList<>(Arrays.asList(cbColumns));
			needed.add("sublibrary");
			for (String col : needed)
				if (!h.hasField(col)) throw new IllegalArgumentException("Column "+col+" not found in "+g.cellsFile+" (expected a called cells table with columns sublibrary,"+String.join(",", cbColumns)+")");
			for (String[] a : cit.loop()) {
				int s = Integer.parseInt(a[h.get("sublibrary")].trim());
				for (String col : cbColumns) {
					String cb = a[h.get(col)];
					unionWithUnderscores.put(cb, "");
					String clean = cb.replace("_", "");
					if (cbLength==-1) cbLength = clean.length();
					else if (cbLength!=clean.length()) throw new RuntimeException("Cell barcodes have different lengths: "+cb);
					cellsPerSublib.computeIfAbsent(s, x->new HashSet<>()).add(clean);
					g.prefixedCalled.add(sublibraryPrefix(s)+clean);
				}
			}
			if (cellsPerSublib.isEmpty()) throw new RuntimeException("No cells found in "+g.cellsFile);

			@SuppressWarnings("unchecked")
			HashSet<String>[] wl = new HashSet[args.length];
			g.whitelist = wl;
			for (int f=0; f<args.length; f++) {
				wl[f] = cellsPerSublib.get(sublib[f]);
				if (wl[f]==null) wl[f] = new HashSet<>();
			}
			for (Integer s : cellsPerSublib.keySet())
				if (!seenSublib.contains(s))
					Gedi.getLog().warning(g.cellsFile+": called cells for sublibrary "+s+" present, but no bam file given!");

			// barcodes.tsv: take library name from there, and check that it uses prefixed barcodes
			g.name = name;
			String bcFile = FileUtils.getExtensionSibling(g.out, ".barcodes.tsv");
			if (new File(bcFile).exists()) {
				HeaderLine bh = new HeaderLine();
				HashSet<String> libs = new HashSet<>();
				long total = 0, found = 0;
				for (String[] a : EI.lines(bcFile).header(bh).split("\t").loop()) {
					libs.add(a[bh.get("Library")]);
					total++;
					if (g.prefixedCalled.contains(a[bh.get("Barcode")])) found++;
				}
				Gedi.getLog().info(String.format("%s: %d/%d barcodes are (prefixed) called barcodes of %s", bcFile, found, total, g.cellsFile));
				if (found<total)
					Gedi.getLog().warning(bcFile+" contains barcodes that will not be in the CIT. Are they prefixed with the sublibrary (AA,AC,...)?");
				if (g.name==null && libs.size()==1) {
					g.name = libs.iterator().next();
					Gedi.getLog().info("Using library name "+g.name+" from "+bcFile);
				}
			} else
				Gedi.getLog().warning(bcFile+" not found; create it with make_atlas_tables.R (with sublibrary prefixes)!");
			if (g.name==null)
				g.name = FileUtils.getNameWithoutExtension(g.out);

			Gedi.getLog().info(String.format("Output %s (%s): %d called barcodes from %s", g.out, g.name, g.prefixedCalled.size(), g.cellsFile));
		}

		BamGenomicRegionStorage storage = new BamGenomicRegionStorage(args);
		if (anti)
			storage.setStrandness(Strandness.Antisense);
		storage.setIgnoreVariations(!var);
		storage.setOnlyPrimary(nosec);
		storage.setKeepReadNames(keepIds);
		// same as 10x: barcode=CB+UB, requires UB!="-", maq 255, prefilter for CB in union of all called cells
		storage.set10x(unionWithUnderscores);
		if (minmaq>=0)
			storage.setMinimalAlignmentQuality(minmaq);

		int numCond = storage.getRandomRecord().getNumConditions();
		if (numCond!=args.length)
			throw new RuntimeException("Expected one condition per bam file ("+args.length+"), got "+numCond+" (XR-count in bam header?)");

		// one writer thread per output
		for (Group g : groups) {
			g.storage = new CenteredDiskIntervalTreeStorage<>(g.out, BarcodedAlignedReadsData.class);
			CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData> st = g.storage;
			g.sink = new IterateIntoSink<>(it->st.fill(it));
		}

		int[] barcodeLength = {-1};
		ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> it = storage.ei();
		if (head>0) it = it.head(head);
		if (progress) it = it.progress(new ConsoleProgress(System.err),-1,r->r.toLocationString());

		while (it.hasNext()) {
			ImmutableReferenceGenomicRegion<AlignedReadsData> r = it.next();
			BarcodedAlignedReadsData ard = (BarcodedAlignedReadsData) r.getData();
			for (Group g : groups) {
				BarcodedAlignedReadsData c = collapse(ard, prefixes, g.whitelist, cbLength, barcodeLength, g.kept, g.droppedNotCalled, g.droppedLength);
				if (c!=null)
					g.sink.put(new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), c));
			}
		}

		for (Group g : groups) {
			g.sink.finish();
			g.storage.setMetaData(DynamicObject.from("conditions",
					DynamicObject.from(new Object[] {
							DynamicObject.from("name",g.name)
					})
					));
		}

		for (Group g : groups) {
			System.out.println("Output\t"+g.out);
			System.out.println("Library\t"+g.name);
			System.out.println("Kept reads\t"+g.kept.get());
			System.out.println("Dropped reads (barcode not called in this sublibrary)\t"+g.droppedNotCalled.get());
			System.out.println("Dropped reads (unexpected barcode length)\t"+g.droppedLength.get());
			System.out.println();
		}
	}

	/**
	 * Collapses the per-sublibrary conditions into a single condition, prefixes the barcodes and
	 * removes barcodes that are not called cells of the respective sublibrary.
	 * @return null if nothing remains
	 */
	private static BarcodedAlignedReadsData collapse(BarcodedAlignedReadsData ard, String[] prefixes, HashSet<String>[] whitelist,
			int cbLength, int[] barcodeLength,
			AtomicLong kept, AtomicLong droppedNotCalled, AtomicLong droppedLength) {

		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
		fac.start();
		ArrayList<DnaSequence> bcs = new ArrayList<>();
		for (int d=0; d<ard.getDistinctSequences(); d++) {
			bcs.clear();
			for (int c=0; c<ard.getNumConditions(); c++) {
				for (DnaSequence bc : ard.getBarcodes(d, c)) {
					String s = bc.toString();
					if (s.length()<=cbLength) {
						droppedLength.incrementAndGet();
						continue;
					}
					if (!whitelist[c].contains(s.substring(0, cbLength))) {
						droppedNotCalled.incrementAndGet();
						continue;
					}
					String p = prefixes[c]+s;
					if (barcodeLength[0]==-1) barcodeLength[0] = p.length();
					if (p.length()!=barcodeLength[0]) {
						droppedLength.incrementAndGet();
						continue;
					}
					bcs.add(new DnaSequence(p));
				}
			}
			if (bcs.size()>0) {
				fac.add(ard, d, v->v, false, false);
				fac.setCount(fac.getDistinctSequences()-1, 0, bcs.size(), bcs.toArray(new DnaSequence[0]));
				kept.addAndGet(bcs.size());
			}
		}
		return fac.getDistinctSequences()==0?null:fac.createBarcode();
	}

	private static void usage() {
		System.out.println("ParseSeqBam2CIT [-p] [-id] [-novar] [-nosec] [-anti] [-minmaq <MAQ>] [-head <n>] [-name <library>] [-cbColumns CB_hex,CB_dT] [-sublibraryPattern <regex>]\n"
				+ "     -cells <called_cells.tsv> <output.cit> <bam1> <bam2> ...\n"
				+ "  or -split <a.called_cells.tsv,b.called_cells.tsv,...> <bam1> <bam2> ...\n\n"
				+ " Converts the bam files of all sublibraries of a Parse Biosciences experiment (STARsolo, CB/UB tags) into a CIT with one condition.\n"
				+ " The cell barcodes are prefixed by two nucleotides identifying the sublibrary (1->AA, 2->AC, 3->AG, 4->AT, 5->CA, ...).\n"
				+ " Only barcodes of called cells (of the respective sublibrary) are kept.\n\n"
				+ " -p shows progress\n"
				+ " -id add ids to CIT\n"
				+ " -cells tsv with columns sublibrary and the cell barcode columns (e.g. auto_cell_calling/called_cells_barcodes.tsv)\n"
				+ " -split comma-separated called cells tables X.called_cells.tsv; writes X.cit for each in one pass over the bams\n"
				+ " -cbColumns comma-separated columns of the cells table containing barcodes (default: CB_hex,CB_dT)\n"
				+ " -sublibraryPattern regex with one group to infer the sublibrary number from the bam path (default: sublibrary_(\\d+))\n"
				+ " -name library/condition name (default: unique Library in <output>.barcodes.tsv, otherwise output file name)\n"
				+ " -minmaq minimal mapping quality (default 255, i.e. unique mappers)\n");
	}

}
