package gedi.util.io.text.fasta;

import java.util.function.Function;

import gedi.util.StringUtils;

public class FastaEntry {
	
	private String sequence;
	private String header;
	private int width = 80;

	public FastaEntry() {
	}
	
	public FastaEntry(String header, String sequence) {
		setHeader(header);
		setSequence(sequence);
	}

	public String getSequence() {
		return sequence;
	}

	public void setSequence(String sequence) {
		this.sequence = sequence;
	}

	/**
	 * Includes the >
	 * @return
	 */
	public String getHeader() {
		return header;
	}

	public void setHeader(String header) {
		this.header = header.startsWith(">")?header:">"+header;
	}
	
	public void setWidth(int width) {
		this.width = width;
	}
	
	public int getWidth() {
		return width;
	}
	
	@Override
	public String toString() {
		if (width<=0)
			return header+"\n"+sequence;
		else
			return header+"\n"+StringUtils.concat("\n", StringUtils.split(sequence,width));
	}
	
	public <T> void headerFromArray(T[] a, int start, int length, String separator, Function<T,String> stringer) {
		StringBuilder sb = new StringBuilder();
		sb.append(">");
		for (int i=start; i<start+length; i++) {
			sb.append(stringer.apply(a[i]));
			sb.append(separator);
		}
		if (sb.length()>1)
			sb.delete(sb.length()-separator.length(), sb.length());
		setHeader(sb.toString());
	}

	public static class ToHeaderTransformer implements Function<FastaEntry,String> {
		@Override
		public String apply(FastaEntry fe) {
			return fe.getHeader();
		}
	}
	
	public static class ToSequenceTransformer implements Function<FastaEntry,String> {
		@Override
		public String apply(FastaEntry fe) {
			return fe.getSequence();
		}
	}
	
	public FastaEntry clone() {
		return new FastaEntry(header,sequence);
	}
	
}
