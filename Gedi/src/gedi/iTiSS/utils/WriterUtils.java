package gedi.iTiSS.utils;

import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

import java.io.IOException;

public class WriterUtils {
    public static LineWriter[] createWriters(String prefix, String[] fileNames, String postfix) {
        LineWriter[] writers = new LineWriter[fileNames.length];
        for (int i = 0; i < fileNames.length; i++) {
            writers[i] = new LineOrientedFile(prefix + fileNames[i] + postfix).write();
        }
        return writers;
    }

    public static void writeLineToAll(LineWriter[] writers, String line) throws IOException {
        for (LineWriter writer : writers) {
            writer.writeLine(line);
        }
    }

    public static void closeWriters(LineWriter[] writers) {
        EI.wrap(writers).forEachRemaining(w -> {
            try {
                w.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
    }
}
