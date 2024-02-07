package gedi.riboseq.analysis;

import java.io.IOException;
import java.util.List;

import gedi.util.io.text.LineWriter;

public interface PriceAnalysis<R> {

	void header(LineWriter out) throws IOException;
	R createContext();
	void reduce(List<R> ctx, LineWriter out);
	
	
	default void plot(String data, String prefix) throws IOException {}
	
	void process(MajorIsoform data, LineWriter out, R ctx);

}
