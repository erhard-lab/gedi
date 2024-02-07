package gedi.util.userInteraction.progress;

import java.util.Iterator;
import java.util.LinkedHashMap;

import gedi.util.functions.EI;
import gedi.util.math.stat.RandomNumbers;

public class ProgressManager implements Runnable {

	private enum ProgressState {
		BeforeFirstView,Running,AfterFinished
	}
	
	private LinkedHashMap<Progress,ProgressState> activeProgresses = new LinkedHashMap<>();
	private Thread thread;
	private boolean started;
	private boolean finishing = false;
	
	private synchronized Thread getThread() {
		synchronized (this) {
			if (thread!=null && started)
				return thread;
			thread = new Thread(this);
			thread.setName("ProgressManager");
			thread.setDaemon(true);
			return thread;	
		}
		
	}

	public void started(Progress p) {
		while (finishing)
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
			}
		
		synchronized (this) {
			activeProgresses.put(p,ProgressState.BeforeFirstView);
		}
		getThread();
		if (!started) {
			thread.start();
			started = true;
		}
	}
	
	public void finished(Progress p) {
		synchronized (this) {
			if (activeProgresses.get(p)==ProgressState.BeforeFirstView)
				activeProgresses.remove(p);
			else
				activeProgresses.put(p,ProgressState.AfterFinished);
			
			finishing = EI.wrap(activeProgresses.values()).filter(s->s==ProgressState.AfterFinished).count()==activeProgresses.size();
			
		}
		
		if (finishing)
			try {
				if (thread!=null)
					thread.join();
				started = false;
				thread = null;
				
			} catch (InterruptedException e) {
			}
		
		synchronized (this) {
			finishing=false;
		}
	}

	@Override
	public void run() {
		
		for (;;) {
			
			synchronized (this) {

				if (activeProgresses.size()==0) break;
				
				for (Progress p : EI.wrap(activeProgresses.keySet()).filter(p->activeProgresses.get(p)==ProgressState.BeforeFirstView).loop()) {
					p.firstView(activeProgresses.size());
					activeProgresses.put(p, ProgressState.Running);
				}
				int index = 0;
				for (Progress p : EI.wrap(activeProgresses.keySet()).loop()) {
					p.updateView(index, activeProgresses.size());
					index++;
				}
				
				Iterator<Progress> it = activeProgresses.keySet().iterator();
				while (it.hasNext()) {
					Progress p = it.next();
					if (activeProgresses.get(p)==ProgressState.AfterFinished) {
						p.lastView(activeProgresses.size());
						it.remove();
					}
				}
				
				activeProgresses.values().removeIf(s->s==ProgressState.AfterFinished);
				
			}
			
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
				return;
			}
		}
		
	}
	
	
	public static void main(String[] args) throws InterruptedException {
		
		System.out.println("\u2639");
		
		ConsoleProgress pr = new ConsoleProgress();
		pr.init();
		pr.setDescription("Singleton");
		
		pr.setCount(100);
		for (int s=0; s<100; s++) {
			Thread.sleep(50);
			pr.incrementProgress();
		}
		
		pr.finish();

		RandomNumbers rnd = RandomNumbers.getGlobal();
		
		ProgressManager man = new ProgressManager();
		for (int i=0;i<5; i++) {
			int ui = i;
			new Thread(()->{
				
				try {
					Thread.sleep((long) (rnd.getUnif()*4000));
					
					ConsoleProgress p = new ConsoleProgress(man);
					p.init();
					p.setDescription("Thread "+ui);
					long pause = (long) (rnd.getUnif()*100);
					
					p.setCount(100);
					for (int s=0; s<100; s++) {
						Thread.sleep(pause);
						p.incrementProgress();
					}
					
					p.finish();
					
				} catch (Exception e) {
				}
				
				
				
			}).start();
		}
	}
	
}
