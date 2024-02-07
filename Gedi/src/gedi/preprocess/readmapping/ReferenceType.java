package gedi.preprocess.readmapping;

public enum ReferenceType {

	Genomic(3,false),Transcriptomic(2,true),Both(4,false),rRNA(1,true),contaminant(1,false),spike(3,false);
	
	public int prio;
	public boolean norc;
	private ReferenceType(int prio, boolean norc) {
		this.prio = prio;
		this.norc = norc;
	}
	
}
