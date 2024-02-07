package gedi.util.program.parametertypes;


import java.util.Map;

import gedi.util.ParseUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

public class ChoicesParameterType<E> implements GediParameterType<E> {

	private Map<String,E> choices;
	private Class<E> cls;
	
	public ChoicesParameterType(Map<String,E> choices, Class<E> cls) {
		this.choices = choices;
		this.cls = cls;
	}

	@Override
	public E parse(String s) {
		E re = ParseUtils.parseChoicesByPrefix(s, true, choices);
		if (re==null) throw new RuntimeException("Could not parse "+s+" (Available: "+EI.wrap(choices.keySet()).concat(",")+")");
		return re;
	}

	@Override
	public Class<E> getType() {
		return cls;
	}
	
	@Override
	public E getDefaultValue() {
		return choices.get(choices.keySet().iterator().next());
	}
	
	@Override
	public boolean parsesMulti() {
		return false;
	}
	
	@Override
	public String helpText() {
		return "Choices are : "+EI.wrap(choices.keySet()).concat(",");
	}

}
