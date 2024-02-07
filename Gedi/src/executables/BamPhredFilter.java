package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;

import gedi.app.Gedi;
import gedi.util.FileUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.r.RRunner;
import gedi.util.sequence.MismatchString;
import gedi.util.sequence.MismatchString.Mismatch;
import gedi.util.sequence.MismatchString.MismatchStringPart;
import gedi.util.sequence.MismatchString.StretchOfMatches;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;


public class BamPhredFilter {

	
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
		if (args.length<1) {
			usage();
			System.exit(1);
		}
		
		Gedi.startup(false);
		
		
		int phred = 33;
		int min = 28;
		Progress progress = new NoProgress();
		String in = null;
		String out = null;
		boolean stats = false;

		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = new ConsoleProgress(System.err);
			else if (args[i].equals("-phred"))
				phred = checkIntParam(args, ++i);
			else if (args[i].equals("-min"))
				min = checkIntParam(args,++i);
			else if (args[i].equals("-stats"))
				stats = true;
			else if (args[i].startsWith("-"))
				throw new IllegalArgumentException("Parameter "+args[i]+" unknown!");
			else {
				in = args[i++];
				if (!new File(in).exists()) throw new RuntimeException("Input file does not exist!");
				if (i<args.length) out = args[i];
				if (i>args.length) {
					usage();
					System.exit(1);
				}
			}
		}
			
		boolean overwrite = out==null;
		if (overwrite) out = FileUtils.getExtensionSibling(in, "tmp");

		SamReader sam = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(in));
		SAMFileWriter output = new SAMFileWriterFactory().makeBAMWriter(sam.getFileHeader(), false, new File(out));
		progress.init();
		
		int maxpos = 0;
		long[][][][] statArr = stats?new long[2][5][5][1024]:null;
		SAMRecordIterator it = sam.iterator();
		while (it.hasNext()) {
			SAMRecord rec = it.next();
			
			MismatchStringPart[] mm = new MismatchString(rec.getStringAttribute("MD")).toArray(new MismatchStringPart[0]);
			int pos = 0;
			int len = mm.length;
			for (int i=0; i<len; i++) {
				if (mm[i] instanceof StretchOfMatches) {
					StretchOfMatches m = (StretchOfMatches) mm[i];
					pos+=m.getGenomicLength();
				} else if (mm[i] instanceof Mismatch) {
					Mismatch m = (Mismatch)mm[i];
					int score = rec.getBaseQualityString().charAt(pos)-phred;
					int filtered = 0;
					if (score<min) {
						if (i==0 || !(mm[i-1] instanceof StretchOfMatches) || i==mm.length-1|| !(mm[i+1] instanceof StretchOfMatches))
							throw new RuntimeException("Malformed mismatch string: "+rec.getStringAttribute("MD"));
						pos+=mm[i+1].getGenomicLength();
						// merge the 3 elements: create new stretch of matches, copy potential subsequent elements, adapt length and decrement i
						mm[i-1] = new StretchOfMatches(mm[i-1].getGenomicLength()+1+mm[i+1].getGenomicLength());
						if (i+2<mm.length) 
							System.arraycopy(mm, i+2, mm, i, mm.length-i-2);
						len-=2;
						i--;
						filtered = 1;
					}
					if (stats) {
						char genomic = m.getGenomicChar();
						char read = rec.getReadString().charAt(pos);
						statArr[filtered][SequenceUtils.inv_nucleotides[genomic]][SequenceUtils.inv_nucleotides[read]][pos]++;
					}
					
					pos+=1;
				}
			}
			maxpos = Math.max(pos, maxpos);
			
			if (len!=mm.length) {
				StringBuilder sb = new StringBuilder();
				for (int i=0; i<len; i++)
					sb.append(mm[i].toString());
				rec.setAttribute("MD", sb.toString());
				int mms = (mm.length-len)/2;
				Integer nm = rec.getIntegerAttribute("nM");
				if (nm!=null) rec.setAttribute("nM", nm-mms);
				
				nm = rec.getIntegerAttribute("NM");
				if (nm!=null) rec.setAttribute("NM", nm-mms);
				
				rec.setAttribute("MD", sb.toString());
				
			}
			output.addAlignment(rec);
			progress.incrementProgress().setDescription(()->rec.getReadName()+" "+rec.getReferenceName()+":"+rec.getAlignmentStart());
		}

		progress.finish();
		output.close();
		sam.close();

		if (overwrite) {
			new File(out).renameTo(new File(in));
			out=in;
		}
		
		if (stats) {
			String tsv = FileUtils.getExtensionSibling(out, "stats.tsv");
			LineWriter st = new LineOrientedFile(tsv).write();
			st.writeLine("Position\tGenomic\tRead\tOriginal\tRetained");
			
			for (int pos=0; pos<maxpos; pos++) {
				for (int g=0; g<4; g++) for (int r=0; r<4; r++) if (r!=g) {
					long ret = statArr[0][g][r][pos];
					long ori = ret + statArr[1][g][r][pos];
					st.writef("%d\t%s\t%s\t%d\t%d\n", pos,SequenceUtils.nucleotides[g],SequenceUtils.nucleotides[r],ori,ret);
				}
			}
			st.close();
			
			try {
				String prefix = FileUtils.getNameWithoutExtension(out);
				File outputFolder = new File(prefix+".stats.plots");
				outputFolder.mkdirs();
				
				Gedi.getLog().info("Running R script for plotting");
				RRunner r2 = new RRunner(out+".phred.stats.R");
				r2.set("prefix",prefix);
				r2.addSource(BamPhredFilter.class.getResourceAsStream("/resources/R/phred.stats.R"));
				r2.run(true);
				
				String json = "{\"plots\":[{\"section\":\"Phred-filter\",\"id\":\"%s_phredfilter\",\"title\":\"%s\","+
								"\"description\":\"Retained mismatches after phred score filtering\",\"img\":\"%s.phred.A>C.png\","+
								"\"imgs\":%s,\"script\":null,\"csv\":\"%s\"}]}";
				StringBuilder imgs = new StringBuilder();
				imgs.append("[");
				for (int g=0; g<4; g++) for (int r=0; r<4; r++) if (r!=g) {
					imgs.append(String.format("{\"name\":\"%s>%s\",\"img\":\"%s.phred.%s>%s.png\"},",SequenceUtils.nucleotides[g],SequenceUtils.nucleotides[r],prefix,SequenceUtils.nucleotides[g],SequenceUtils.nucleotides[r]));
				}
				imgs.deleteCharAt(imgs.length()-1);
				imgs.append("]");
				
				json = String.format(json, prefix.replaceAll("\\.", "_"),prefix,prefix,imgs.toString(),tsv);
				FileUtils.writeAllText(json, new File(outputFolder,prefix+".phred.report.json"));
			} catch (Throwable e) {
				Gedi.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		
	}


	private static void usage() {
		System.out.println("BamPhredFilter [-stats] [-p] [-phred <offset>] [-min <score>] <in.bam> [<out.bam>]\n\n -stats generate stats file\n -p shows progress\n -phred <offset> ascii to score conversion offset, either 33 or 64 (default 33)\n -min <score> (minimal phred score to keep mismatch (default: 28)\n\nIf no out.bam is given, the input is overwritten!");
	}
	
}
