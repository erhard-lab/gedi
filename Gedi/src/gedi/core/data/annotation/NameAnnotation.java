package gedi.core.data.annotation;


import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class NameAnnotation implements BinarySerializable, NameProvider, Comparable<NameAnnotation> {
	
	private String name;
	
	public NameAnnotation() {
	}
	
	public NameAnnotation(String name) {
		this.name = name;
	}

	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putString(name);
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		name = in.getString();
	}
	
	@Override
	public int compareTo(NameAnnotation o) {
		return getName().compareTo(o.getName());
	}
	
	@Override
	public String toString() {
		return name;
	}

	public String getName() {
		return name;
	}
	
	public void setName(String name) {
		this.name = name;
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
		NameAnnotation other = (NameAnnotation) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		return true;
	}

	
	
}
