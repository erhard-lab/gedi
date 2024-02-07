package gedi.util.program;

import java.io.File;
import java.util.ArrayList;

import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GediParameterType;

public class GediParameter<T> {

	private T defaultValue;
	private T value;
	private ArrayList<T> values;
	private GediParameterType<T> type;
	private String name;
	private String description;
	private GediParameterSet set;
	boolean optional;
	boolean multi;
	boolean removeFile = false;
	private String[] shortcut = null; 

	public GediParameter(GediParameterSet set, String name, String description, boolean multi, GediParameterType<T> type) {
		this(set, name,description,multi,type,type.getDefaultValue(),false);
	}
	
	public GediParameter(GediParameterSet set, String name, String description, boolean multi, GediParameterType<T> type, T defaultValue) {
		this(set,name,description,multi,type,defaultValue,false);
	}
	
	public GediParameter(GediParameterSet set, String name, String description, boolean multi, GediParameterType<T> type, boolean optional) {
		this(set, name,description,multi,type,type.getDefaultValue(),optional);
	}
	
	public GediParameter(GediParameterSet set, String name, String description, boolean multi, GediParameterType<T> type, T defaultValue, boolean optional) {
		this.set = set;
		this.value = this.defaultValue = defaultValue;
		this.name = name;
		this.description = description;
		this.multi = multi;
		this.type = type;
		this.optional = optional;
		if (multi) values = new ArrayList<>();
		set.add(this);
		
		if (type instanceof FileParameterType)
			setRemoveFile(((FileParameterType) type).isTemp());
	}

	public GediParameter<T> setShortcut(String... shortcut) {
		this.shortcut = shortcut;
		return this;
	}
	
	public boolean isShortcut() {
		return shortcut!=null;
	}
	
	public String[] getShortcut() {
		return shortcut;
	}
	
	public String getName() {
		return name;
	}
	
	public GediParameter<T> setRemoveFile(boolean removeFile) {
		this.removeFile = removeFile;
		return this;
	}
	
	public boolean isRemoveFile() {
		return removeFile;
	}
	
	public boolean isOptional() {
		return optional;
	}
	
	public boolean isMulti() {
		return multi;
	}
	
	public String getUsage(int pad) {
		StringBuilder sb = new StringBuilder();
		if (isOptional()) sb.append("[");
		
		if (isFile()) sb.append(getFile());
		else {
			sb.append(" -").append(name);
			if (type.hasValue())
				sb.append(" <").append(name).append(">");
		}
		
		if (pad==-1) return sb.toString();
		
		pad-=sb.length();
		for (int i=0; i<pad; i++)
			sb.append(' ');
		sb.append(description);
		
		if (defaultValue!=null && type.hasValue()) sb.append(" (default: ").append(defaultValue).append(")");
		
		String help = type.helpText();
		if (help!=null)
			sb.append(" ").append(help);
		
		if (isOptional()) sb.append("]");
		
		return sb.toString();
	}
	
	public boolean hasDefault(){
		return defaultValue!=null;
	}
	
	public String getDescription() {
		return description;
	}
	
	public void set(T value) {
		if (isMulti())
			this.values.add(value);
		this.value = value;
	}
	
	public boolean isFile() {
		return type.getType()==File.class;
	}
	
	public File getFile() {
		return FileParameterType.getFile(getName(),getParameterSet());
	}

	public T get() {
		if (isFile()) {
			File f = getFile();
			if (f.exists() && (f.length()>0 || (f.isDirectory() && EI.fileNames(f.getPath()).count()>0))) return (T) f;
			return null;
		}
		return value;
	}
	
	public String getStringDescriptor() {
		if (isMulti())
			return EI.wrap(values).concat(",");
		return StringUtils.toString(value);
	}
	
	public T getValue() {
		return value;
	}
	
	public ArrayList<T> getList() {
		if (values.isEmpty() && defaultValue!=null) values.add(defaultValue);
		return values;
	}
	
	public GediParameterType<T> getType() {
		return type;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GediParameter other = (GediParameter) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	@Override
	public String toString() {
		return name+"="+get();
	}

	public GediParameterSet getParameterSet() {
		return set;
	}

	

	
	
	
	
}
