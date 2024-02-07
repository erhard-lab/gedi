package gedi.grand3.processing;

import java.io.File;
import java.io.IOException;
import java.util.function.IntFunction;

import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.util.ArrayUtils;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public class SubreadStatistics implements SubreadCounter<SubreadStatistics> {

	@Override
	public SubreadStatistics spawn(int index) {
		return new SubreadStatistics(subreads.length,reindex);
	}

	@Override
	public void integrate(SubreadStatistics other) {
		for (int c=0; c<lengths.length; c++)
				ArrayUtils.add(lengths[c], other.lengths[c]);
		for (int c=0; c<subreads.length; c++)
			ArrayUtils.add(subreads[c], other.subreads[c]);
	}
	
	
	// cat,subread,mismatch,condition
	private double[][] lengths;
	private double[][] subreads;
	private int[] reindex;
	
	public SubreadStatistics(int numSubreads, int[] reindex) {
		lengths = new double[1001][ArrayUtils.max(reindex)+1];
		subreads = new double[numSubreads][ArrayUtils.max(reindex)+1];
		this.reindex = reindex;
	}

	@Override
	public void count(SubreadProcessorMismatchBuffer buffer) {
		int len = Math.min(buffer.getRead().getRegion().getTotalLength(),lengths.length-1);
		buffer.count(lengths[len], reindex);
		
		SubreadsAlignedReadsData sr = buffer.getRead().getData();
		int[] lens = new int[subreads.length];
		for (int i=0; i<sr.getNumSubreads(buffer.distinct); i++) 
			lens[sr.getSubreadId(buffer.distinct, i)]+=sr.getSubreadEnd(buffer.distinct, i, buffer.getRead().getRegion().getTotalLength())-sr.getSubreadStart(buffer.distinct, i);
		
		for (int s=0; s<subreads.length; s++) {
			buffer.count(subreads[s], lens[s], reindex);
		}
		
	}
	
	public void writeLengths(File out, IntFunction<String> condNameFunction) throws IOException {
		LineWriter wr = new LineOrientedFile(out.getPath()).write();
		wr.writeLine("Condition\tLength\tFrequency");
		for (int cond=0; cond<lengths[0].length; cond++)
			for (int l=0; l<lengths.length; l++)
				if (lengths[l][cond]>0)
					wr.writef("%s\t%d\t%.0f\n", 
							condNameFunction.apply(cond),
							l,
							lengths[l][cond]
									);
		
		wr.close();
	}
	

	public void writeSubreads(File out, IntFunction<String> condNameFunction, String[] semantic) throws IOException {
		LineWriter wr = new LineOrientedFile(out.getPath()).write();
		wr.writeLine("Condition\tSubread\tFrequency");
		for (int cond=0; cond<lengths[0].length; cond++)
			for (int l=0; l<subreads.length; l++)
				if (subreads[l][cond]>0)
					wr.writef("%s\t%s\t%.0f\n", 
							condNameFunction.apply(cond),
							semantic[l],
							subreads[l][cond]
									);
		
		wr.close();
	}
	

	
}

