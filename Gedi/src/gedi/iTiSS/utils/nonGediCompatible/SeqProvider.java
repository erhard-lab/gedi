package gedi.iTiSS.utils.nonGediCompatible;

import gedi.core.region.GenomicRegion;
import gedi.core.sequence.SequenceProvider;

import java.util.HashSet;
import java.util.Set;

public class SeqProvider implements SequenceProvider {
    private int seqLength;
    private Set<String> seqNames;

    public SeqProvider(String seqName, int seqLength) {
        this.seqLength = seqLength;
        seqNames = new HashSet<>();
        seqNames.add(seqName+"+");
        seqNames.add(seqName+"-");
        seqNames.add(seqName);
    }

    @Override
    public Set<String> getSequenceNames() {
        return seqNames;
    }

    @Override
    public CharSequence getPlusSequence(String name, GenomicRegion region) {
        return null;
    }

    @Override
    public char getPlusSequence(String name, int pos) {
        return 0;
    }

    @Override
    public int getLength(String name) {
        return seqNames.contains(name) ? seqLength : -1;
    }
}
