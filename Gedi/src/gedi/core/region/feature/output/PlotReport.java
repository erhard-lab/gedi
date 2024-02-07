package gedi.core.region.feature.output;

public class PlotReport {
	public String section;
	public String id;
	public String title;
	public String description;
	public String img;
	public SecondaryPlot[] imgs;
	public String script;
	public String csv;
	public PlotReport(String section,String id, String title, String description, String img, SecondaryPlot[] imgs, String script, String csv) {
		this.section = section;
		this.id = id;
		this.title = title;
		this.description = description;
		this.img = img;
		this.imgs = imgs;
		this.script = script;
		this.csv = csv;
	}
	
	@Override
	public String toString() {
		return title;
	}
	
}