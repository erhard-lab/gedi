package gedi.util.datastructure.array;

import gedi.util.functions.IntDoubleConsumer;
import gedi.util.functions.IntDoubleToDoubleFunction;

public interface IndexDoubleProcessor {

	
	/**
	 * Iterates over all non-zero entries (ordered), applies op, and changes the contents!
	 * @param op
	 */
	default void process(IntDoubleToDoubleFunction op) {
		process(op,0,length());
	}
	void process(IntDoubleToDoubleFunction op, int start, int end);
	
	/**
	 * Iterates over all non-zero entries (ordered) and applies op.
	 * @param op
	 */
	default void iterate(IntDoubleConsumer op) {
		iterate(op,0,length());
	}
	default void iterate(IntDoubleConsumer op, int start, int end) {
		process((pos,val)->{
			op.accept(pos, val);
			return val;
		},start,end);
	}

	int length();

	default double[] toArray() {
		double[] re = new double[length()];
		iterate((index,value)->re[index]=value);
		return re;
	}
	default double[] toArray(int start, int end) {
		double[] re = new double[end-start];
		iterate((index,value)->re[index-start]=value,start,end);
		return re;
	}
}
