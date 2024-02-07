package gedi.util.job.schedule;

import gedi.util.job.ExecutionContext;

import java.util.concurrent.ExecutorService;
import java.util.function.Consumer;

public interface PetriNetScheduler extends Runnable {

	ExecutionContext getExecutionContext();
	void addListener(PetriNetListener rg);
	void removeListener(PetriNetListener rg);
	void setFinishAction(Consumer<ExecutionContext> finishAction);
	
	
}
