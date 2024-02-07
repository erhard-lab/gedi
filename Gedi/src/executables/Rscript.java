package executables;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.logging.Logger;

import gedi.app.Gedi;
import gedi.app.classpath.ClassPathCache;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.mutable.MutablePair;
import gedi.util.r.RRunner;

public class Rscript {

	private static final Logger log = Logger.getLogger( Rscript.class.getName() );
	public static void main(String[] args) {
		try {
			start(args);
		} catch (UsageException e) {
			usage("An error occurred: "+e.getMessage(),null);
			if (ArrayUtils.find(args, "-D")>=0)
				e.printStackTrace();
		} catch (Exception e) {
			System.err.println("An error occurred: "+e.getMessage());
			if (ArrayUtils.find(args, "-D")>=0)
				e.printStackTrace();
		}
	}
	
	private static void usage(String message, String additional) {
		System.err.println();
		if (message!=null){
			System.err.println(message);
			System.err.println();
		}
		System.err.println("R <Options> script.R...");
		System.err.println();
		System.err.println("Options:");
		System.err.println(" --x.y[3].z=val\t\tDefine a template variable; val can be json (e.g. --x='{\\\"prop\\\":\\\"val\\\"}'");
		System.err.println(" -j <json-file>\t\tDefine template variables");
		System.err.println(" -r <r-file>\t\tKeep the script file with this path (implies -merge)");
		System.err.println(" -merge\t\tMerge all discovered scripts into one file");
		System.err.println(" -dontrun \t\tDo not run");
		System.err.println();
		System.err.println(" -D\t\t\tOutput debugging information");
		System.err.println(" -h\t\t\tShow this message");
		System.err.println();
		if (additional!=null) {
			System.err.println(additional);
			System.err.println("");
		}
	}
	
	
	private static class UsageException extends Exception {
		public UsageException(String msg) {
			super(msg);
		}
	}
	
	
	private static String checkParam(String[] args, int index) throws UsageException {
		if (index>=args.length || args[index].startsWith("-")) throw new UsageException("Missing argument for "+args[index-1]);
		return args[index];
	}
	
	private static String[] checkPair(String[] args, int index) throws UsageException {
		int p = args[index].indexOf('=');
		if (!args[index].startsWith("--") || p==-1) throw new UsageException("Not an assignment parameter (--name=val): "+args[index]);
		return new String[] {args[index].substring(2, p),args[index].substring(p+1)};
	}
	
	
	
	public static void start(String[] args) throws Exception {
		Gedi.startup(true);
		
		DynamicObject param = DynamicObject.getEmpty();
		String script = null;
		boolean dontrun = false;
		boolean merge = false;
		
		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-h")) {
				usage(null,null);
				return;
			}
			else if (args[i].equals("-D")) {
			}else if (args[i].equals("-j")) {
				String path = checkParam(args, ++i);
				DynamicObject param1 = DynamicObject.parseJson(FileUtils.readAllText(new File(path)));
				param = DynamicObject.cascade(param,param1);
			}else if (args[i].equals("-r")) {
				script = checkParam(args, ++i);
				merge = true;
			}else if (args[i].equals("-dontrun")) {
				dontrun = true;
			}else if (args[i].equals("-merge")) {
				merge = true;
			}
			else if (args[i].startsWith("--")) {
				String[] p = checkPair(args,i);
				
				DynamicObject param1=DynamicObject.parseExpression(p[0], DynamicObject.parseJsonOrString(p[1]));
				param = DynamicObject.cascade(param,param1);
			}
			else if (!args[i].startsWith("-")) 
					break;
			else throw new UsageException("Unknown parameter: "+args[i]);
			
		}

		
		boolean deleteScript = script==null;
		ArrayList<MutablePair<String, InputStream>> rscripts = new ArrayList<>();
		
		for (; i<args.length; i++) {
			
			// try to find rscript:
			String rscript = args[i];
			if (Rscript.class.getResourceAsStream(rscript)!=null) {
				rscripts.add(new MutablePair<>(rscript,Rscript.class.getResourceAsStream(rscript)));
				log.info("Running "+rscript);
			}
			else if (Rscript.class.getResourceAsStream("/resources/R/"+rscript)!=null) {
				rscripts.add(new MutablePair<>("/resources/R/"+rscript,Rscript.class.getResourceAsStream("/resources/R/"+rscript)));
			}
			else if (Rscript.class.getResourceAsStream("/resources/"+rscript)!=null) {
				rscripts.add(new MutablePair<>("/resources/"+rscript,Rscript.class.getResourceAsStream("/resources/"+rscript)));
			}
			else if (new File(rscript).exists() && new File(rscript).isFile()) {
				rscripts.add(new MutablePair<>(rscript,new FileInputStream(new File(rscript))));
			}
			else {
				int found = 0;
				for (String rr : ClassPathCache.getInstance().getResourcesByExtension("R").str().matches(rscript).loop()) {
					rscripts.add(new MutablePair<>(rr,ClassPathCache.getInstance().getClassPathOfFile(rr).getResourceAsStream(rr)));
					found++;
				}
				if (found==0)
					throw new UsageException("Script not found!");
			}
		}
		
		
		
		if (merge) {
			if (script==null) script = File.createTempFile("Rscript", ".R").getAbsolutePath();
			RRunner r = new RRunner(script);
			for (String n : param.getProperties())
				r.set(n,param.getEntry(n));
			for (MutablePair<String,InputStream> p : rscripts) {
				log.info("Merging "+p.Item1);
				r.addSource(p.Item2);
			}
			log.info("Running ...");
			if (dontrun)
				r.dontrun(deleteScript);
			else
				r.run(deleteScript);
			if (!deleteScript) log.info("Script file: "+script);
		} else {
			for (MutablePair<String,InputStream> p : rscripts) {
				String tscript = script;
				if (tscript!=null && rscripts.size()>1) tscript = FileUtils.getExtensionSibling(tscript, FileUtils.getNameWithoutExtension(p.Item1)+"."+FileUtils.getExtension(tscript));
				if (tscript==null) tscript = File.createTempFile("Rscript", ".R").getAbsolutePath();
				RRunner r = new RRunner(tscript);
				for (String n : param.getProperties())
					r.set(n,param.getEntry(n));

				log.info("Running "+p.Item1+"...");
				r.addSource(p.Item2);
				if (dontrun)
					r.dontrun(deleteScript);
				else
					r.run(deleteScript);
				if (!deleteScript) log.info("Script file: "+tscript);
			}
		}

	}

	
}
