package gedi.util.functions;

import java.util.concurrent.LinkedBlockingQueue;
import java.util.function.Consumer;

/**
 * Put is thread-safe!
 * @author erhard
 *
 * @param <T>
 */
public class IterateIntoSink<T> implements Consumer<T> {

	private LinkedBlockingQueue<T> queue;
	private T terminator = (T) new Object();
	private Thread t;
	
	public IterateIntoSink(Consumer<ExtendedIterator<T>> sink) {
		this(sink,1024);
	}
	public IterateIntoSink(Consumer<ExtendedIterator<T>> sink, int buffer) {
		this.queue = new LinkedBlockingQueue<T>(buffer);
		
		
		t = new Thread() {
			@Override
			public void run() {
				try {
					sink.accept(new SinkIterator());
				} catch (Exception e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
		};
		t.setDaemon(true);
		t.setName("SinkIterator");
		t.start();
		
	}
	
	
	public void put(T element) throws InterruptedException {
		queue.put(element);
	}
	
	public void finish() throws InterruptedException {
		queue.put(terminator);
		t.join();
	}

	
	private class SinkIterator implements ExtendedIterator<T> {
		private T n;
	
		@Override
		public boolean hasNext() {
			if(n==null)
				try {
					n=queue.take();
				} catch (InterruptedException e) {
					return false;
				}
			return n!=terminator;
		}
	
	
		@Override
		public T next() {
			if(n==null)
				try {
					n=queue.take();
				} catch (InterruptedException e) {
				}
			T re = n;
			n = null;
			return re;
		}
	}


	@Override
	public void accept(T t) {
		try {
			put(t);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
	}
	
	
}
