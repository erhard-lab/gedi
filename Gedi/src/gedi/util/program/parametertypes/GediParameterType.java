package gedi.util.program.parametertypes;

import java.io.File;

import gedi.util.program.GediParameterSet;

public interface GediParameterType<T> {

	
	T parse(String s);
	Class<T> getType();
	default boolean hasValue() {
		return true;
	}
	default T getDefaultValue() {
		return null;
	}
	default boolean isInternal() {
		return false;
	}
	default boolean parsesMulti() {
		return false;
	}
	default String helpText() {
		return null;
	}
	
}
