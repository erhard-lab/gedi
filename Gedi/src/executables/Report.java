package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;

import javax.script.ScriptException;

import gedi.core.region.feature.output.PlotReport;
import gedi.util.FileUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.io.text.LineIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.StreamLineWriter;
import gedi.util.io.text.jhp.Jhp;

public class Report {
	
	private String[] files;
	
	public Report(String[] files) {
		this.files = files;
	}
	
	
	public void write(LineWriter writer) throws IOException, ScriptException {
		DynamicObject d = DynamicObject.getEmpty();
		for (String l : files)
			d = d.merge(DynamicObject.parseJson(FileUtils.readAllText(new File(l))));
		
		PlotReport[] plots = d.getEntry("plots").javafy(PlotReport[].class);

		LinkedHashMap<String, ArrayList<PlotReport>> sections = new LinkedHashMap<>();
		for (PlotReport p : plots) 
			sections.computeIfAbsent(p.section, x->new ArrayList<>()).add(p);
		
		
		Jhp jhp = new Jhp();
		jhp.getJs().setSelf(d);
		if (!d.hasProperty("title"))
			jhp.getJs().putVariable("title", "Report");
		jhp.getJs().putVariable("sections", sections);
		
		String html = new LineIterator(Report.class.getResource("/resources/templates/stats/plots.html.jhp").openStream()).concat("\n");
		html = jhp.apply(html);
		
		writer.write(html);
		writer.close();
	}

	public static void main(String[] args) throws IOException, ScriptException {
		if ((args.length==1 && args[0].equals("-h")) || args.length==0) {
			System.err.println("Usage: Report <json> [<json> ...]");
			System.exit(1);
		}
		
		new Report(args).write(new StreamLineWriter(System.out));
	}
	
}
