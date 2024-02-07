package gedi.core.sequence;

import gedi.core.region.GenomicRegion;
import gedi.util.functions.EI;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class CompositeSequenceProvider implements SequenceProvider {

	private ArrayList<SequenceProvider> providers = new ArrayList<SequenceProvider>();

	public CompositeSequenceProvider add(SequenceProvider provider) {
		providers.add(provider);
		return this;
	}
	
	@Override
	public CharSequence getPlusSequence(String name, GenomicRegion region) {
		for (SequenceProvider seq :providers){
			CharSequence re = seq.getPlusSequence(name, region);
			if (re!=null) return re;
		}
		return null;
	}

	
	@Override
	public char getPlusSequence(String name, int pos) {
		for (SequenceProvider seq :providers){
			char re = seq.getPlusSequence(name, pos);
			if (re!=0) return re;
		}
		return '\0';
	}


	@Override
	public Set<String> getSequenceNames() {
		return EI.wrap(providers).demultiplex(s->s.getSequenceNames().iterator()).toCollection(new HashSet<>());
	}

	@Override
	public int getLength(String name) {
		for (SequenceProvider seq :providers){
			int re = seq.getLength(name);
			if (re!=-1) return re;
		}
		return -1;
	}
	
	public ArrayList<SequenceProvider> getProviders() {
		return providers;
	}
	
}
