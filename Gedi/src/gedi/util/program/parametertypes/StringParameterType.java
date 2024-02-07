package gedi.util.program.parametertypes;

import gedi.util.parsing.StringParser;

public class StringParameterType implements GediParameterType<String> {

	@Override
	public String parse(String s) {
		return new StringParser().apply(s);
	}

	@Override
	public Class<String> getType() {
		return String.class;
	}

}
