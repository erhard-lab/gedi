package gedi.util.program.parametertypes;


import gedi.core.genomic.Genomic;
import gedi.util.functions.EI;

public class GenomicParameterType implements GediParameterType<Genomic> {

	@Override
	public Genomic parse(String s) {
		return Genomic.get(EI.split(s, ' '));
	}

	@Override
	public Class<Genomic> getType() {
		return Genomic.class;
	}
	
	@Override
	public boolean parsesMulti() {
		return true;
	}

}
