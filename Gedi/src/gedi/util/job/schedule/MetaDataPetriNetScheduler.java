package gedi.util.job.schedule;

import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.collections.longcollections.LongArrayList;
import gedi.util.dynamic.DynamicObject;
import gedi.util.job.ExecutionContext;
import gedi.util.job.FireTransition;
import gedi.util.job.Place;
import gedi.util.job.Transition;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.mutable.MutableLong;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.locks.LockSupport;
import java.util.function.Consumer;
import java.util.logging.Level;
import java.util.logging.Logger;

public class MetaDataPetriNetScheduler implements PetriNetScheduler {

	private static final Logger log = Logger.getLogger( PetriNetScheduler.class.getName() );
	
	protected ExecutionContext context;
	protected Consumer<ExecutionContext> finishAction;
	
	private boolean logging = false;
	private ArrayList<PetriNetListener> listeners = new ArrayList<PetriNetListener>();
	
	public MetaDataPetriNetScheduler(ExecutionContext context) {
		this.context = context;
	}
	
	
	@Override
	public void addListener(PetriNetListener rg) {
		listeners.add(rg);
	}
	
	@Override
	public void removeListener(PetriNetListener rg) {
		listeners.remove(rg);
	}
	
	@Override
	public void setFinishAction(Consumer<ExecutionContext> finishAction) {
		this.finishAction = finishAction;
	}
	
	public void setLogging(boolean logging) {
		this.logging = logging;
	}

	@Override
	public void run() {
		
		try {
			int eid = context.startExecution();
			
			PetriNetEvent event = new PetriNetEvent(eid, context);
			
			for (PetriNetListener l : listeners)
				l.petriNetExecutionStarted(event);
			
			long start = System.nanoTime();
			
			
			HashSet<Transition> ready = new HashSet<Transition>();
			for (Place p : context.getPetrNet().getSources())
				addReadyConsumers(p, ready);
			ArrayList<Transition> iter = new ArrayList<Transition>();
			
			HashSet<Future<FireTransition>> futures = new HashSet<Future<FireTransition>>();
			if (Thread.interrupted() || context.getExecutionId()!=eid) {
				if (logging) log.log(Level.FINE,"Canceled execution "+eid);
				for (PetriNetListener l : listeners)
					l.petriNetExecutionCancelled(event);
				return;
			}
						
			if (logging) log.log(Level.FINE,"Start executing Petri net (id="+eid+") Disabled: "+context.getDisabledTransitions());
			
			while (!ready.isEmpty()) {
				
				iter.clear();
				iter.addAll(ready);
				ready.clear();
				
				for (Transition n : iter) {
					if (Thread.interrupted() || context.getExecutionId()!=eid) {
						for (Future<FireTransition> f : futures)
							f.cancel(true);
						if (logging) log.log(Level.INFO,"Canceled execution "+eid);
						for (PetriNetListener l : listeners)
							l.petriNetExecutionCancelled(event);
						return;
					}
					
					if (logging) log.log(Level.FINE,"Submitting "+n+" (id="+eid+") "+context);
					DynamicObject meta = n.getJob().meta(context.createMeta(n));
					context.putToken(n.getOutput(), meta);
					if (logging) log.log(Level.FINE,"Finished "+n+" (id="+eid+")");
					addReadyConsumers(n.getOutput(), ready);
					
				}
			}
			
			
			if (logging) log.log(Level.FINE,String.format("Finished executing Petri net (id="+eid+") in %d ns\n",System.nanoTime()-start));
			
			for (PetriNetListener l : listeners)
				l.petriNetExecutionFinished(event);
			
			if (finishAction!=null)
				finishAction.accept(context);
			
		} catch (Throwable e) {
			log.log(Level.SEVERE,"Uncaught exception!",e);
			
		}
	}

	private void addReadyConsumers(Place p, HashSet<Transition> ready) {
		for (Transition t : p.getConsumers())
			if (!ready.contains(t) && context.isReady(t))
				ready.add(t);
	}


	@Override
	public ExecutionContext getExecutionContext() {
		return context;
	}

	
	
	
}
