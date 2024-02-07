package gedi.util.nashorn;

import java.util.function.Function;

import javax.script.ScriptException;

import gedi.util.mutable.MutableTuple;
import org.openjdk.nashorn.api.scripting.ScriptObjectMirror;

public class JSFunction<I,O> implements Function<I,O> {

	private ScriptObjectMirror p;
	private boolean useAsThis;
	
	public JSFunction(String code) throws ScriptException {
		this(false,code);
	}
	public JSFunction(boolean useAsThis, String code) throws ScriptException {
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
	
	
	@SuppressWarnings("unchecked")
	@Override
	public O apply(I t) {
		if (t instanceof MutableTuple)
			return (O) p.call(null, ((MutableTuple)t).getArray());
		if (useAsThis)
			return (O) p.call(t);
		return (O) p.call(null, t);
	}

}
