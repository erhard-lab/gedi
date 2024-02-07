package executables;

import java.util.ArrayList;
import java.util.logging.Logger;

import gedi.app.Gedi;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.io.text.jhp.TemplateEngine;

public class Template {

	private static final Logger log = Logger.getLogger( Template.class.getName() );
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
		System.err.println("Template <Options>");
		System.err.println();
		System.err.println("Options:");
		System.err.println(" --x.y[3].z=val\t\tDefine a template variable; val can be json (e.g. --x='{\\\"prop\\\":\\\"val\\\"}'");
		System.err.println(" -j [json-file]\t\tDefine template variables");
		System.err.println(" -t [template-file]\t\tProcess a template");
		System.err.println("");
		System.err.println(" -h [<template>]	Print usage information (of this program and the given template)");
		System.err.println();
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
	
	
	private static int checkMultiParam(String[] args, int index, ArrayList<String> re) throws UsageException {
		while (index<args.length && !args[index].startsWith("-")) 
			re.add(args[index++]);
		return index-1;
	}
	private static String checkParam(String[] args, int index) throws UsageException {
		if (index>=args.length || args[index].startsWith("-")) throw new UsageException("Missing argument for "+args[index-1]);
		return args[index];
	}
	
	private static int checkIntParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new UsageException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(args[index]);
	}

	private static double checkDoubleParam(String[] args, int index) throws UsageException {
		String re = checkParam(args, index);
		if (!StringUtils.isNumeric(re)) throw new UsageException("Must be a double: "+args[index-1]);
		return Double.parseDouble(args[index]);
	}
	private static String[] checkPair(String[] args, int index) throws UsageException {
		int p = args[index].indexOf('=');
		if (!args[index].startsWith("--") || p==-1) throw new UsageException("Not an assignment parameter (--name=val): "+args[index]);
		return new String[] {args[index].substring(2, p),args[index].substring(p+1)};
	}
	
	public static void start(String[] args) throws Exception {
		Gedi.startup(true);
		
		TemplateEngine te = new TemplateEngine();

		int i;
		for (i=0; i<args.length; i++) {
			
			if (args[i].equals("-h")) {
				ArrayList<String> list = new ArrayList<>();
				i = checkMultiParam(args, i+1, list);
				if (list.isEmpty())
					usage(null,null);
				else {
					usage(null,"\n"+te.help(list.get(0)));
				}
				return;
			}
			else if (args[i].equals("-D")) {
			}else if (args[i].equals("-j")) {
				te.json(checkParam(args, ++i));
			}
			else if (args[i].startsWith("--")) {
				String[] p = checkPair(args,i);
				te.parameter(p[0], p[1]);
			}else if (args[i].equals("-t")) {
				te.template(checkParam(args, ++i));
			}
			else throw new UsageException("Unknown parameter: "+args[i]);
			
		}

		
		System.out.println(te.toString());
		
	}
	
}
