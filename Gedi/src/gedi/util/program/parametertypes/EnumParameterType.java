package gedi.util.program.parametertypes;


import gedi.util.ParseUtils;
import gedi.util.StringUtils;

public class EnumParameterType<E extends Enum<E>> implements GediParameterType<E> {

	private Class<E> enumClass;
	
	public EnumParameterType(Class<E> enumClass) {
		this.enumClass = enumClass;
	}

	@Override
	public E parse(String s) {
		E re = ParseUtils.parseEnumNameByPrefix(s, true, enumClass);
		if (re==null) throw new RuntimeException("Could not parse "+s+" as a "+enumClass.getName()+" (Available: "+ParseUtils.getEnumTrie(enumClass, true).keySet().toString()+")");
		return re;
	}

	@Override
	public Class<E> getType() {
		return enumClass;
	}
	
	@Override
	public E getDefaultValue() {
		return enumClass.getEnumConstants()[0];
	}
	
	@Override
	public boolean parsesMulti() {
		return false;
	}
	
	@Override
	public String helpText() {
		return "Choices are : "+StringUtils.concat(",", enumClass.getEnumConstants());
	}

}
