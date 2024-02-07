package gedi.util.nashorn;

import java.util.function.Predicate;

import javax.script.ScriptException;

import gedi.util.mutable.MutableTuple;
import org.openjdk.nashorn.api.scripting.ScriptObjectMirror;

public class JSPredicate<T> implements Predicate<T> {


	private ScriptObjectMirror p;
	private boolean useAsThis;
	
	public JSPredicate(String code) throws ScriptException {
		this(false,code);
	}
	public JSPredicate(boolean useAsThis, String code) throws ScriptException {
		this.useAsThis = useAsThis;
		
		if (!code.startsWith("function(")) {
			StringBuilder c = new StringBuilder();
			c.append("function() {\n");
			if (code.contains(";")) {
				c.append(code);
				c.append("}");
			} else {
				c.append("return "+code+";\n}");
			}
			code = c.toString();
		}
		p = new JS().execSource(code);
	}
	
	
	@Override
	public boolean test(T t) {
		if (useAsThis)
			return (boolean) p.call(t);
		if (t instanceof MutableTuple)
			return (boolean) p.call(null, ((MutableTuple)t).getArray());
		return (boolean) p.call(null, t);
	}

	
}
