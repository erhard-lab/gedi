package gedi.grand3.experiment;

import gedi.util.SequenceUtils;

public class MetabolicLabel {
	
	public enum MetabolicLabelType {
		_4sU {
			@Override
			public char getGenomic() {
				return 'T';
			}
			@Override
			public char getRead() {
				return 'C';
			}
		},
		_6sG {
			@Override
			public char getGenomic() {
				return 'G';
			}
			@Override
			public char getRead() {
				return 'A';
			}
		};
		
		public abstract char getGenomic();
		public abstract char getRead();
		public char getGenomicReverse() {
			return SequenceUtils.getDnaComplement(getGenomic());
		}
		public char getReadReverse() {
			return SequenceUtils.getDnaComplement(getRead());
		}

		public String toString() {
			return name().substring(1);
		}
		
		public static MetabolicLabelType fromString(String s) {
			return MetabolicLabelType.valueOf("_"+s);
		}
	}
	

	private MetabolicLabelType type;
	private double concentration;
	private double duration;
	private double chase;
	public MetabolicLabel(MetabolicLabelType type, double concentration, double duration, double chase) {
		super();
		this.type = type;
		this.concentration = concentration;
		this.duration = duration;
		this.chase = chase;
	}
	public MetabolicLabelType getType() {
		return type;
	}
	public double getConcentration() {
		return concentration;
	}
	public double getDuration() {
		return duration;
	}
	public double getChase() {
		return chase;
	}
	@Override
	public String toString() {
		return "MetabolicLabel [type=" + type + ", concentration=" + concentration + ", duration=" + duration
				+ ", chase=" + chase + "]";
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		long temp;
		temp = Double.doubleToLongBits(chase);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(concentration);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		temp = Double.doubleToLongBits(duration);
		result = prime * result + (int) (temp ^ (temp >>> 32));
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
		MetabolicLabel other = (MetabolicLabel) obj;
		if (Double.doubleToLongBits(chase) != Double.doubleToLongBits(other.chase))
			return false;
		if (Double.doubleToLongBits(concentration) != Double.doubleToLongBits(other.concentration))
			return false;
		if (Double.doubleToLongBits(duration) != Double.doubleToLongBits(other.duration))
			return false;
		if (type != other.type)
			return false;
		return true;
	}

	
	
}