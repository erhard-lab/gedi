package gedi.util.io.text;


public class StringLineWriter implements LineWriter, CharSequence {

	
	private StringBuilder sb = new StringBuilder();

	
	@Override
	public void write(String line) {
		sb.append(line);
	}

	@Override
	public void flush() {
	}

	@Override
	public void close() {
	}

	
	@Override
	public String toString() {
		return sb.toString();
	}

	public void setContent(String content) {
		sb.delete(0, sb.length());
		sb.append(content);
	}

	@Override
	public int length() {
		return sb.length();
	}

	@Override
	public char charAt(int index) {
		return sb.charAt(index);
	}

	@Override
	public CharSequence subSequence(int start, int end) {
		return sb.subSequence(start, end);
	}
	
}
