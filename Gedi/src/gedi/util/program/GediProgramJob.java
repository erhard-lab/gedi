package gedi.util.program;


import java.util.logging.Logger;

import gedi.util.job.ExecutionContext;
import gedi.util.mutable.MutableTuple;

public class GediProgramJob extends DummyJob {

	private GediProgram s;
	private Class[] input;
	private boolean dry;
	
	public GediProgramJob(GediProgram s) {
		this.s = s;
		input = new Class[s.getInputSpec().size()];
		for (int i=0; i<input.length; i++)
			input[i] = Boolean.class;
	}
	
	public GediProgramJob setDry(boolean dry) {
		this.dry = dry;
		return this;
	}

	@Override
	public Class[] getInputClasses() {
		return input;
	}
	
	@Override
	public Boolean execute(ExecutionContext context, MutableTuple input) {
		if (dry) {
			Logger log = ((GediProgramContext)context.getContext("context")).getLog();
			synchronized (log) {
				log.info("Running "+s.getName()+"  In: "+s.getInputSpec().getNames()+" Out: "+s.getOutputSpec().getNames());	
			}
			return true;
		}
		
		try {
			GediProgramContext pc = context.getContext("context");
			pc.startMessage(s.getName());
			s.execute(pc);
			pc.finishMessage(s.getName());
		} catch (Exception e) {
			throw new RuntimeException("Could not run "+s,e);
		}
		return true;
	}

	@Override
	public String getId() {
		return s.getName()+ind;
	}
	
	private boolean disabled = false;
	
	@Override
	public boolean isDisabled(ExecutionContext context) {
		return disabled;
	}
	
	public void disable() {
		this.disabled = true;
	}

}
