package gedi.util.program.parametertypes;

import gedi.util.parsing.DoubleParser;

public class DoubleParameterType implements GediParameterType<Double> {

	@Override
	public Double parse(String s) {
		return new DoubleParser().apply(s);
	}

	@Override
	public Class<Double> getType() {
		return Double.class;
	}

}
