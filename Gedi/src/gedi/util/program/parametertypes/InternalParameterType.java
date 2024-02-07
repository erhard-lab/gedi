package gedi.util.program.parametertypes;


public class InternalParameterType<T> implements GediParameterType<T> {

	private Class<T> cls;
	public InternalParameterType(Class<T> cls) {
		this.cls = cls;
	}

	@Override
	public T parse(String s) {
		throw new RuntimeException("Not possible!");
	}

	@Override
	public Class<T> getType() {
		return cls;
	}

	@Override
	public boolean isInternal() {
		return true;
	}
	
}
