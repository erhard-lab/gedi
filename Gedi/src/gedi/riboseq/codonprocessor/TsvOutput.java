package gedi.riboseq.codonprocessor;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableTuple;

public class TsvOutput extends CodonProcessorOutput {

	@Override
	public void createOutput2(String[] conditions) throws IOException {
		LineWriter out = new LineOrientedFile(counter.getPrefix()+".tsv").write();
		for (int i=0; i<counter.getVariableNames().length; i++) {
			out.write(counter.getVariableNames()[i]);
			out.write("\t");
		}
		out.writeLine(EI.wrap(conditions).concat("\t"));
		
		
		ArrayList<MutableTuple> keys = new ArrayList<MutableTuple>(counter.keys());
		Collections.sort(keys);
		
		for (MutableTuple key : keys) {
			NumericArray n = counter.getCounts(key);
			for (int i=0; i<counter.getVariableNames().length; i++) {
				out.write(key.get(i).toString());
				out.write("\t");
			}
				
			for (int c=0; c<conditions.length; c++) {
				if (c>0) out.write("\t");
				out.write(n.format(c));
			}
			
			out.writeLine();
		}
		
		out.close();
	}

	
}
