package gedi.util.functions;

import java.util.function.Consumer;

@FunctionalInterface
public interface Consumer2<T> extends Consumer<T>{

	
	 void accept2(T t) throws Exception;
	 
	 default void accept(T t) {
		 try {
			accept2(t);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	 }
}
