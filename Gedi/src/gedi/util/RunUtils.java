package gedi.util;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.ProcessBuilder.Redirect;
import java.util.function.Consumer;

public class RunUtils {
	
	
	public static String output(String...cmds)  {
		try {
			File tmp = File.createTempFile("run", ".out");
			Process proc = new ProcessBuilder(cmds)
				.redirectError(Redirect.INHERIT)
				.redirectOutput(Redirect.to(tmp))
				.redirectInput(Redirect.INHERIT)
				.start();
			proc.waitFor();
			
			String re = FileUtils.readAllText(tmp);
			tmp.delete();
			return re;
		} catch (IOException | InterruptedException e) {
			throw new RuntimeException("Cannot run "+StringUtils.concat(" ", cmds),e);
		}
	}

	
	public static int run(String...cmds)  {
		try {
			Process proc = new ProcessBuilder(cmds)
				.redirectError(Redirect.INHERIT)
				.redirectOutput(Redirect.INHERIT)
				.redirectInput(Redirect.INHERIT)
				.start();
			return proc.waitFor();
		} catch (IOException | InterruptedException e) {
			throw new RuntimeException("Cannot run "+StringUtils.concat(" ", cmds),e);
		}
	}
	
	public static int pipeInto(Consumer<PrintWriter> input, String...cmds)  {
		try {
			Process proc = new ProcessBuilder(cmds)
				.redirectError(Redirect.INHERIT)
				.redirectOutput(Redirect.INHERIT)
				.start();
			PrintWriter wr = new PrintWriter(proc.getOutputStream());
			input.accept(wr);
			wr.close();
			return proc.waitFor();
		} catch (IOException | InterruptedException e) {
			throw new RuntimeException("Cannot run "+StringUtils.concat(" ", cmds),e);
		}
	}

}
