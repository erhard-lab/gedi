package gedi.iTiSS.utils;

import gedi.core.genomic.Genomic;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.functions.EI;
import org.apache.commons.math3.distribution.BinomialDistribution;

import java.util.*;

public class SequenceUtils2 {
    public static final Set<String> startCodons = new HashSet<String>();
    static {
        startCodons.add("ATG");
        startCodons.add("CTG");
        startCodons.add("GTG");
        startCodons.add("ACG");
        startCodons.add("ATC");
        startCodons.add("TTG");
        startCodons.add("ATT");
        startCodons.add("ATA");
        startCodons.add("AGG");
        startCodons.add("AAG");
    }

    private static List<StringBuilder> duplicateList(List<StringBuilder> lst) {
        List<StringBuilder> newlst = new ArrayList<>();
        for (int l = 0; l < lst.size(); l++) {
            StringBuilder newsb = new StringBuilder();
            newsb.append(lst.get(l).toString());
            newlst.add(newsb);
        }
        return newlst;
    }

    private static void addToLst(List<StringBuilder> lst, char ch) {
        for (int w = 0; w < lst.size(); w++) {
            lst.get(w).append(ch);
        }
    }

    public static boolean hasMotif(ReferenceGenomicRegion<?> inRgr, Genomic genomic, String seq, int min, int max) {
        return getMotifPosition(inRgr, genomic, seq, min, max) != Integer.MAX_VALUE;
    }

    public static List<String> motifToSequences(String motif) {
        List<StringBuilder> seqsToTest = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        seqsToTest.add(sb);
        for (int s = 0; s < motif.length(); s++) {
            if (motif.charAt(s)=='W') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'T');
                addToLst(dup, 'A');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='Y') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'C');
                addToLst(dup, 'T');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='R') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'A');
                addToLst(dup, 'G');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='B') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                List<StringBuilder> dup2 = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'C');
                addToLst(dup, 'G');
                addToLst(dup2, 'T');
                seqsToTest.addAll(dup);
                seqsToTest.addAll(dup2);
            } else if (motif.charAt(s)=='S') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'C');
                addToLst(dup, 'G');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='K') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'T');
                addToLst(dup, 'G');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='M') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'A');
                addToLst(dup, 'C');
                seqsToTest.addAll(dup);
            } else if (motif.charAt(s)=='D') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                List<StringBuilder> dup2 = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'A');
                addToLst(dup, 'G');
                addToLst(dup2, 'T');
                seqsToTest.addAll(dup);
                seqsToTest.addAll(dup2);
            } else if (motif.charAt(s)=='V') {
                List<StringBuilder> dup = duplicateList(seqsToTest);
                List<StringBuilder> dup2 = duplicateList(seqsToTest);
                addToLst(seqsToTest, 'A');
                addToLst(dup, 'C');
                addToLst(dup2, 'G');
                seqsToTest.addAll(dup);
                seqsToTest.addAll(dup2);
            } else {
                addToLst(seqsToTest, motif.charAt(s));
            }
        }
        return EI.wrap(seqsToTest).map(StringBuilder::toString).list();
    }

    public static double getMotifExpProb(String motif, int min, int max) {
        List<String> seqs = motifToSequences(motif);
        double prob = seqs.size()/Math.pow(4.0, motif.length());
        int searchSpace = (max+1)-min;
        int maxOccurrences = searchSpace-(motif.length()-1);
        BinomialDistribution binom = new BinomialDistribution(maxOccurrences, prob);
        double cumProb = 0;
        for (int i = 1; i <= maxOccurrences; i++) {
            cumProb += binom.probability(i);
        }
        return cumProb;
    }

    public static int getMotifPosition(ReferenceGenomicRegion<?> inRgr, Genomic genomic, String seq, int min, int max) {
        List<String> seqsToTest = motifToSequences(seq);
        if (max > 0) {
            inRgr = new ImmutableReferenceGenomicRegion<>(inRgr.getReference(), inRgr.getRegion().union(inRgr.getDownstream(max).getRegion()));
        }
        inRgr = new ImmutableReferenceGenomicRegion<>(inRgr.getReference(), inRgr.getRegion().union(inRgr.getUpstream(min*-1).getRegion()));
        for (String actualSeq : seqsToTest) {
            int seqIndex = genomic.getSequence(inRgr).toString().substring(0,(max+1)-min).indexOf(actualSeq);
            if (seqIndex >= 0) {
//                return (max - (actualSeq.toString().length()-seqIndex))+1;
                return min+seqIndex;
            }
        }
        return Integer.MAX_VALUE;
    }

    public static boolean isStartCodon(String codon, boolean allowAlternativeStartCodons) {
        if (!allowAlternativeStartCodons) {
            return codon.equals("ATG");
        }
        return startCodons.contains(codon);
    }

    public static boolean codonEqualsAny(String codon, String... codonsToCheck) {
        for (String c : codonsToCheck) {
            if (codon.equals(c)) {
                return true;
            }
        }
        return false;
    }
}
