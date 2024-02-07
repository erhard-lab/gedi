package gedi.proteomics.maxquant;

import java.io.IOException;

import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class MassspecEvidence implements BinarySerializable {

	private String type;
	private int evidenceCount;
	private int multiplicity;
	public MassspecEvidence(String type, int evidenceCount, int multiplicity) {
		super();
		this.type = type;
		this.evidenceCount = evidenceCount;
		this.multiplicity = multiplicity;
	}
	public String getType() {
		return type;
	}
	public int getEvidenceCount() {
		return evidenceCount;
	}
	public int getMultiplicity() {
		return multiplicity;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + evidenceCount;
		result = prime * result + multiplicity;
		result = prime * result + ((type == null) ? 0 : type.hashCode());
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
		MassspecEvidence other = (MassspecEvidence) obj;
		if (evidenceCount != other.evidenceCount)
			return false;
		if (multiplicity != other.multiplicity)
			return false;
		if (type == null) {
			if (other.type != null)
				return false;
		} else if (!type.equals(other.type))
			return false;
		return true;
	}
	@Override
	public String toString() {
		return type + "[n=" + evidenceCount + ", a=" + multiplicity
				+ "]";
	}
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		out.putString(type);
		out.putCInt(evidenceCount);
		out.putCInt(multiplicity);
	}
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		this.type = in.getString();
		this.evidenceCount = in.getCInt();
		this.multiplicity = in.getCInt();
	}
	
	
	
}
