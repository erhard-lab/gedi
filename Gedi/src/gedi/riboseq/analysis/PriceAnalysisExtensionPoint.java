package gedi.riboseq.analysis;

import gedi.app.extension.DefaultExtensionPoint;

public class PriceAnalysisExtensionPoint extends DefaultExtensionPoint<String, PriceAnalysis>{

	
	private static PriceAnalysisExtensionPoint instance;

	public static PriceAnalysisExtensionPoint getInstance() {
		if (instance==null) 
			instance = new PriceAnalysisExtensionPoint(PriceAnalysis.class);
		return instance;
	}

	protected PriceAnalysisExtensionPoint(
			Class<PriceAnalysis> cls) {
		super(cls);
	}


}
