package gedi.riboseq.codonprocessor;

import java.io.IOException;
import java.util.function.Consumer;

public abstract class CodonProcessorOutput {

	protected CodonProcessorCounter counter;
	
	public void setCounter(CodonProcessorCounter counter) {
		if (this.counter!=null) throw new RuntimeException("This output is already attached to another counter!");
		this.counter = counter;
	}
	
	
	public void createOutput(String[] conditions) {
		try {
			createOutput2(conditions);
		} catch (IOException e) {
			throw new RuntimeException("Could not write output for "+counter.getPrefix()+"!",e);
		}
	}
	
	public abstract void createOutput2(String[] conditions) throws IOException;

}
