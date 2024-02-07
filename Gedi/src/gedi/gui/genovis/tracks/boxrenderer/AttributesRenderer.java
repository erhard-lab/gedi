package gedi.gui.genovis.tracks.boxrenderer;

import gedi.core.data.annotation.AttributesProvider;
import gedi.core.data.annotation.NameProvider;
import gedi.util.PaintUtils;
import gedi.util.StringUtils;
import gedi.util.dynamic.DynamicObject;

import java.awt.Color;

public class AttributesRenderer extends BoxRenderer<Object> {

	
	private DynamicObject features;
	private String attribute = "name";

	public AttributesRenderer() {
		stringer = (rgr)->{
			if (features==null) return "";
			String t = getType(rgr.getData());
			if (t==null) return "";
			String s = features.getEntry(t).getEntry("label").asString();
			return StringUtils.replaceVariables(s,p->{
				Object r = rgr.getData() instanceof AttributesProvider?((AttributesProvider)rgr.getData()).getAttribute(p):null;
				return r==null?p:r.toString();
			});
		};
		background = (rgr)->{
			if (features==null) return Color.gray;
			String t = getType(rgr.getData());
			if (t==null) return Color.gray;
			String s = features.getEntry(t).getEntry("color").asString();
			return PaintUtils.parseColor(s);
		};
	}
	
	
	private String getType(Object a) {
		if (attribute.equals("name") && a instanceof NameProvider)
			return ((NameProvider)a).getName();
		if (a instanceof AttributesProvider)
			return ((AttributesProvider) a).getStringAttribute(attribute);
		return null;
	}


	public void setAttribute(String attribute) {
		this.attribute = attribute;
	}
	
	public void setFeatures(DynamicObject features) {
		this.features = features;
	}
	
	


}
