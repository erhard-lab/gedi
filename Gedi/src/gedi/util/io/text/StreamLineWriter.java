package gedi.util.io.text;

import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;

public class StreamLineWriter implements LineWriter {

	
	private Writer writer;

	public StreamLineWriter(OutputStream stream) {
		writer = new OutputStreamWriter(stream);
	}
	
	public StreamLineWriter(Writer writer) {
		this.writer = writer;
	}
	
	
	@Override
	public void write(String line) throws IOException {
		writer.write(line);
	}

	@Override
	public void flush() throws IOException {
		writer.flush();
	}

	@Override
	public void close() throws IOException {
		writer.close();
	}

	
	
}
