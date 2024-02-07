package gedi.util.functions;

@FunctionalInterface
public interface IntObjectConsumer<T> {

	void accept(int n, T o);
	
}
