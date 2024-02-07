package gedi.util.nashorn;

import java.util.function.Consumer;

import javax.script.ScriptException;

import gedi.util.mutable.MutableTuple;
import org.openjdk.nashorn.api.scripting.ScriptObjectMirror;

public class JSConsumer<I> implements Consumer<I> {

	private ScriptObjectMirror p;
	private boolean useAsThis;
	
	public JSConsumer(String code) throws ScriptException {
		this(false,code);
	}
	public JSConsumer(boolean useAsThis, String code) throws ScriptException {
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
	public void accept(I t) {
		if (t instanceof MutableTuple)
			p.call(null, ((MutableTuple)t).getArray());
		else if (useAsThis)
			p.call(t);
		else
			p.call(null, t);
	}

}
