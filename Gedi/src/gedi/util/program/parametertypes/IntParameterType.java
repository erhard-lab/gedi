package gedi.util.program.parametertypes;

import gedi.util.parsing.IntegerParser;

public class IntParameterType implements GediParameterType<Integer> {

	@Override
	public Integer parse(String s) {
		return new IntegerParser().apply(s);
	}

	@Override
	public Class<Integer> getType() {
		return Integer.class;
	}

}
