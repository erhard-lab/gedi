package gedi.util.functions;

import java.util.function.BiConsumer;

@FunctionalInterface
public interface BiConsumer2<T,S> extends BiConsumer<T,S>{

	
	 void accept2(T t,S s) throws Exception;
	 
	 default void accept(T t,S s) {
		 try {
			accept2(t,s);
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	 }
}
