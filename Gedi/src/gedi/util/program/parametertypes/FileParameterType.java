package gedi.util.program.parametertypes;

import java.io.File;

import gedi.util.StringUtils;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;

public class FileParameterType implements GediParameterType<File> {
	
	
	private boolean isTemp = false;
	
	public FileParameterType() {
	}
	
	public FileParameterType(boolean isTemp) {
		this.isTemp = isTemp;
	}

	public boolean isTemp() {
		return isTemp;
	}
	
	@Override
	public File parse(String s) {
		return new File(s);
	}

	@Override
	public Class<File> getType() {
		return File.class;
	}

	public static File getFile(String name, GediParameterSet set) {
		return new File(StringUtils.replaceVariables(name, s->{
			GediParameter<Object> param = set.get(s);
			if (param==null)
				throw new RuntimeException("Parameter "+s+" unknown!");
			if (param.get()==null) return "${"+s+"}";
			return StringUtils.toString(param.get());
		}));
	}
	
}
