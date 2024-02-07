package gedi.util.program.parametertypes;

import gedi.util.parsing.LongParser;

public class LongParameterType implements GediParameterType<Long> {

	@Override
	public Long parse(String s) {
		return new LongParser().apply(s);
	}

	@Override
	public Class<Long> getType() {
		return Long.class;
	}

}
