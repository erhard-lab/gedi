package gedi.iTiSS.utils;

import gedi.iTiSS.merger2.TsrFileEntry;
import gedi.util.functions.EI;

import java.util.Comparator;
import java.util.List;

public class TsrFileUtils {
    // 'Best TSR' in this context means the TSR, which was called by the most distinct datasets (greatest getScore())
    // and, if multiple TSRs with the same totalDelta exist, takes the one with the greatest readcount
    public static TsrFileEntry getBestTsr(List<TsrFileEntry> entries) {
        if (entries.size() == 0) {
            return null;
        }
        int maxScore = EI.wrap(entries).mapToInt(TsrFileEntry::getScore).max();
        return EI.wrap(entries).filter(f -> f.getScore() == maxScore).sort(Comparator.comparingDouble(tsrFileEntry -> -tsrFileEntry.getMaxReadCount())).first();
    }
}
