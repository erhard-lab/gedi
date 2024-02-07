package gedi.util.functions;

public interface ParallelizedState<T extends ParallelizedState<T>> {
	
	T spawn(int index);
	/**
	 * Contract: other might be useless afterwards!
	 * @param other
	 */
	void integrate(T other);

}
