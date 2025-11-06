package gedi.iTiSS.utils.datastructures;

import gedi.util.SequenceUtils;
import gedi.util.io.text.LineWriter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

public class SenseToAntisensePairs {
    List<SenseToAntisenseInfoHolder> entries;

    public SenseToAntisensePairs() {
        entries = new ArrayList<>();
    }

    public void sortEntries(Comparator<SenseToAntisenseInfoHolder> comparator) {
        entries.sort(comparator);
    }

    public void writeOutStats(LineWriter writer) throws IOException {
        writer.writeLine(SenseToAntisenseInfoHolder.getStatsHeader());
        for (SenseToAntisenseInfoHolder holder : entries) {
            writer.writeLine(holder.toString());
        }
    }

    public void writeOutSenseFasta(LineWriter writer) throws IOException {
        for (SenseToAntisenseInfoHolder holder : entries) {
            writer.writeLine(">" + holder.getFastaId());
            writer.writeLine(holder.getSequence());
        }
    }

    public void writeOutAntenseFasta(LineWriter writer) throws IOException {
        for (SenseToAntisenseInfoHolder holder : entries) {
            writer.writeLine(">" + holder.getFastaId());
            String sequence = SequenceUtils.getDnaReverseComplement(holder.getSequence());
            writer.writeLine(sequence);
        }
    }

    public void add(SenseToAntisenseInfoHolder info) {
        entries.add(info);
    }
}
