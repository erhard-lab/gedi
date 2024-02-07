package gedi.app.extension;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;



public abstract class AbstractPrioritiesExtensionPoint<C,K,T extends ToDoubleFunction<C>> implements ExtensionPoint<K,T> {

	
	private Class<T> cls;
	protected HashMap<K,ArrayList<Class<? extends T>>> ext = new HashMap<>();
	private BiFunction<Predicate<K>,C,K> key;
	
	protected AbstractPrioritiesExtensionPoint(Class<T> cls, BiFunction<Predicate<K>,C,K> key) {
		this.cls = cls;
		this.key = key;
	}
	
	@Override
	public Class<T> getExtensionPointClass() {
		return cls;
	}

	@Override
	public void addExtension(Class<? extends T> extension, K key) {
		ext.computeIfAbsent(key, x->new ArrayList<>()).add(extension);
	}

	protected ExtensionContext empty = new ExtensionContext();
	
	public T get(ExtensionContext context, C ck) {
		try {
			if (context==null) context = empty;
			
			ArrayList<Class<? extends T>> l = ext.get(key.apply(ext.keySet()::contains,ck));
			if (l==null) 
				return null;
			
			T re = null;
			double prio = -1;
			for (Class<? extends T> cls : l) {
				T t = context.newInstance(cls);
				double p = t.applyAsDouble(ck);
				if (p>=0 && (re==null || p>prio)) {
					re = t;
					prio = p;
				}
			}
			return re;
		} catch (InstantiationException | IllegalAccessException | IllegalArgumentException | InvocationTargetException e) {
			Logger.getLogger(getClass().getName()).log(Level.SEVERE, "Could not get extension for "+key,e);
			throw new RuntimeException("Could not get extension for "+key,e);
		}
	}
	
	
}
