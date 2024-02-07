package executables;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import gedi.util.FileUtils;
import gedi.util.functions.EI;
import gedi.util.r.RRunner;

public class TestCIT {

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
		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else if (args[i].equals("-c"))
				clear = true;
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
		
		
		File[] readstats = EI.wrap(args).map(c->{
			return new File(new File (new File(c).getParentFile(),"report"),FileUtils.getNameWithoutExtension(c)+".reads.tsv");
		}).toArray(File.class);
		if (EI.wrap(readstats).mapToInt(f->f.exists()?1:0).sum()==args.length) {
			String name = FileUtils.getNameWithoutExtension(out);
			String script = "report/"+name+".reads.R";
			out = new File(new File("report"), name+".reads.tsv").getPath();
			
			RRunner r = new RRunner(script);
			r.setNumeric("files", "c("+EI.wrap(readstats).map(f->"\""+f.getPath()+"\"").concat(",")+")");
			r.set("out", out);
			r.addSource(r.getClass().getResourceAsStream("/resources/R/mergeread.R"));
			r.run(false);
			
			String tsv = out;
			for (String png : EI.fileNames("report").filter(f->f.startsWith(name+".reads") && f.endsWith(".png")).loop()) {
				String title = FileUtils.getNameWithoutExtension(png);
				FileUtils.writeAllText("{\"plots\":[{\"section\":\"Mapping summary\",\"id\":\"ID"+png.replace('.', '_')+"\",\"title\":\""+title+"\",\"description\":\"Mapping summary. For a description see mapping statistics!\",\"img\":\""+png+"\",\"script\":\""+script+"\",\"csv\":\""+tsv+"\"}]}",new File("report/"+FileUtils.getExtensionSibling(png, ".report.json")));
			}

		}
		
		
	}

	private static void usage() {
		System.out.println("TestCIT [-c] [-p] [-s skip1,skip2,...] <output> <file1> <file2> ... \n\n -c removes the input files after successful merging\n -p shows progress\n -s skip chromosomes");
	}
	
}
