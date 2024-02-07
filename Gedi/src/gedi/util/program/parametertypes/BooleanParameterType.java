package gedi.util.program.parametertypes;


public class BooleanParameterType implements GediParameterType<Boolean> {

	@Override
	public Boolean parse(String s) {
		return true;
	}

	@Override
	public Class<Boolean> getType() {
		return Boolean.class;
	}

	@Override
	public boolean hasValue() {
		return false;
	}
	
	@Override
	public Boolean getDefaultValue() {
		return false;
	}
}
