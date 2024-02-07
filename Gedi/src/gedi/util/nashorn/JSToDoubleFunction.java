package gedi.util.nashorn;

import java.util.function.ToDoubleFunction;

import javax.script.ScriptException;

import org.openjdk.nashorn.api.scripting.ScriptObjectMirror;


public class JSToDoubleFunction<I> implements ToDoubleFunction<I> {

	private ScriptObjectMirror p;
	
	public JSToDoubleFunction(String code) throws ScriptException {
		p = new JS().execSource(code);
		
	}
	
	
	@Override
	public double applyAsDouble(I t) {
		return (double) p.call(null, t);
	}

}
