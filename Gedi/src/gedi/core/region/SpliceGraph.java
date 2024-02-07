package gedi.core.region;

import gedi.core.data.annotation.Transcript;
import gedi.util.datastructure.tree.redblacktree.IntervalTreeSet;
import gedi.util.datastructure.tree.redblacktree.SimpleInterval;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

import java.util.function.Consumer;

public class SpliceGraph {

	
	private ArrayGenomicRegion reg;
	
	
	private IntervalTreeSet<Intron> introns = new IntervalTreeSet<Intron>(null);
	
	public SpliceGraph(int start, int end) {
		this.reg = new ArrayGenomicRegion(start, end);
	}
	
	public void addIntron(int start, int end) {
		if (start>reg.getStart() && end<reg.getEnd())
			introns.add(new Intron(start, end-1));
	}
	
	public void addIntrons(GenomicRegion region) {
		for (int i=1; i<region.getNumParts(); i++)
			addIntron(region.getEnd(i-1), region.getStart(i));
	}
	

	public boolean contains(int start, int end) {
		return introns.contains(new Intron(start, end-1));
	}
	
	public ExtendedIterator<Intron> iterateIntrons() {
		return EI.wrap(introns);
	}
	
	public IntervalTreeSet<Intron> getIntrons() {
		return introns;
	}
	

	public void forEachIntronStartingBetween(int start, int end, Consumer<Intron> consumer) {
		introns.forEachIntervalIntersecting(start, end-1, intron->{
			if (intron.getStart()>=start && intron.getStart()<end)
				consumer.accept(intron);
		});
		
	}
	
	
	public void forEachIntronEndingBetween(int start, int end, Consumer<Intron> consumer) {
		introns.forEachIntervalIntersecting(start, end-1, intron->{
			if (intron.getEnd()>=start && intron.getEnd()<end)
				consumer.accept(intron);
		});
		
	}
	
	@Override
	public String toString() {
		return introns.toString();
	}
	
	
	
	
	public static class Intron extends SimpleInterval {

		public Intron(int start, int stop) {
			super(start, stop);
		}

		
	}


}
