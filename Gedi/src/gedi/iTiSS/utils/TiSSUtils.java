package gedi.iTiSS.utils;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTree;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.numeric.GenomicNumericProvider;
import gedi.core.data.numeric.diskrmq.DiskGenomicNumericLoader;
import gedi.core.data.numeric.diskrmq.DiskGenomicNumericProvider;
import gedi.core.data.reads.*;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.*;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.iTiSS.merger.Tiss;
import gedi.iTiSS.merger.Tsr;
import gedi.iTiSS.utils.datastructures.SparseNumericArray;
import gedi.iTiSS.utils.loader.TsrData;
import gedi.iTiSS.utils.machineLearning.PeakAndPos;
import gedi.riboseq.utils.RiboUtils;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;
import java.util.stream.Collectors;

public class TiSSUtils {

    /**
     * Returns a 2D-array containing the indices of each replicate in each row with the possibility
     * to skip certain positions by the {@code skip} character
     * Example:
     * Input: xoxo111222
     * Output: [[0,2], [1,3], [4,5,6], [7,8,9]]
     * @param replicatesString input string
     * @param skip a character indicating positions to skip. \u0000 (the null character) for no skip
     * @return 2D-array with indices of replicates in each row
     */
    public static int[][] extractReplicatesFromString(String replicatesString, char skip) {
        String uniqueReps = StringUtils2.uniqueCharacters(StringUtils.remove(replicatesString, skip));
        int[][] repAry = new int[uniqueReps.length()][];
        for (int i = 0; i < uniqueReps.length(); i++) {
            char current = uniqueReps.charAt(i);
            repAry[i] = StringUtils2.indicesOfChar(replicatesString, current);
        }
        return repAry;
    }

    public static int[][] extractReplicatesFromString(String replicatesString) {
        return extractReplicatesFromString(replicatesString, Character.MIN_VALUE);
    }

    /**
     * Returns a 2D-array containing the indices of each replicate in each row.
     * Replicates are given by {@code reps}.
     * Each row contains the a timepoint indicated by the given characters in {@code timecourseString}.
     * They are ordered by {@code Character}-ordering.
     * @param timecourseString input string
     * @param reps replicate string
     * @param skip a character indicating positions to skip. \u0000 (the null character) for no skip
     * @return 2D-array with indices of replicates in each row
     */
    public static int[][] extractTimecoursesFromString(String timecourseString, int[][] reps, char skip) {
        int[][] tcAry = new int[reps.length][];
        for (int i = 0; i < tcAry.length; i++) {
            String rep = new String(StringUtils2.extract(timecourseString, reps[i]));
            int[] tc = extractSingleTimecourse(rep, skip);
            for (int j = 0; j < tc.length; j++) {
                tc[j] = reps[i][tc[j]];
            }
            tcAry[i] = tc;
        }
        return tcAry;
    }

    private static int[] extractSingleTimecourse(String timecourseString, char skip) {
        if (!hasOnlyUniqueCharacters(StringUtils.remove(timecourseString, skip))) {
            System.out.println("Only unique characters are allowed inside a single replicate for the time points");
            return null;
        }
        String uniqueReps = StringUtils2.uniqueCharacters(StringUtils.remove(timecourseString, skip));
        List<Character> timepoints = new ArrayList<>();
        for (char c : uniqueReps.toCharArray()) {
            timepoints.add(c);
        }
        timepoints.sort(Character::compareTo);
        int[] tcAry = new int[timepoints.size()];

        for (int i = 0; i < tcAry.length; i++) {
            tcAry[i] = StringUtils2.indexOf(timecourseString, timepoints.get(i));
        }

        return tcAry;
    }

    private static boolean hasOnlyUniqueCharacters(String string) {
        return StringUtils2.uniqueCharacters(string).length() == string.length();
    }

    public static int[][] extractTimecoursesFromString(String timecourseString, int[][] reps) {
        return extractTimecoursesFromString(timecourseString, reps, Character.MIN_VALUE);
    }

    public static List<List<Integer>> groupConcurrents(List<Integer> lst, int gap) {
        List<List<Integer>> out = new ArrayList<>();
        if (lst.size() <= 1) {
            out.add(new ArrayList<>(lst));
            return out;
        }

        List<Integer> cpy = new ArrayList<>(lst);
        cpy.sort(Integer::compare);
        List<Integer> concurrents = new ArrayList<>();
        int last = cpy.get(0);
        concurrents.add(last);
        for (int i = 1; i < cpy.size(); i++) {
            if (cpy.get(i) - last <= gap) {
            } else {
                out.add(concurrents);
                concurrents = new ArrayList<>();
            }
            last = cpy.get(i);
            concurrents.add(last);
        }
        out.add(concurrents);
        return out;
    }

    public static Map<ReferenceSequence, List<Integer>> extractTissFromFile(String path, int skip) throws IOException {
        Map<ReferenceSequence, List<Integer>> out = new HashMap<>();
        EI.lines(path).skip(skip).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            ReferenceSequence ref = Chromosome.obtain(split[0]);
            int tiss = Integer.parseInt(split[1]);
            if (!out.containsKey(ref)) {
                out.put(ref, new ArrayList<>());
            }
            out.get(ref).add(tiss);
        });
        return out;
    }

    public static Map<ReferenceSequence, List<List<Integer>>> extractTissFromFileGrouped(String path, int skip, int gap) throws IOException {
        Map<ReferenceSequence, List<Integer>> tiss = extractTissFromFile(path, skip);
        Map<ReferenceSequence, List<List<Integer>>> out = new HashMap<>();
        for (ReferenceSequence k : tiss.keySet()) {
            List<Integer> tissList = tiss.computeIfAbsent(k, t -> new ArrayList<>());
            out.put(k, groupConcurrents(tissList, gap));
        }
        return out;
    }

    public static Map<ReferenceSequence, List<Tsr>> extractTsrsFromFile(String path, int skip, int gap, int id) throws IOException {
        Map<ReferenceSequence, List<Tsr>> tsrMap = new HashMap<>();
        Map<ReferenceSequence, List<List<Integer>>> tissGrouped = extractTissFromFileGrouped(path, skip, gap);
        for(ReferenceSequence ref : tissGrouped.keySet()) {
            List<Tsr> tsrList = EI.wrap(tissGrouped.get(ref)).map(tissGroup -> {
                Tsr tsr = new Tsr(createTiss(id, tissGroup.get(0)));
                for (int i = 1; i < tissGroup.size(); i++) {
                    Set<Integer> calledBy = new HashSet<>();
                    calledBy.add(id);
                    tsr.add(new Tiss(calledBy, tissGroup.get(i)));
                }
                return tsr;
            }).list();
            tsrMap.put(ref, tsrList);
        }
        return tsrMap;
    }

    public static Map<ReferenceSequence, List<Tsr>> extractTsrsFromFinalFile(String path, int skip) throws IOException {
        Map<ReferenceSequence, Map<GenomicRegion, Tsr>> tsrMap = new HashMap<>();
        EI.lines(path).skip(skip).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            Map<GenomicRegion, Tsr> reg2tsr = tsrMap.computeIfAbsent(Chromosome.obtain(split[0]), empty -> new HashMap<>());
            GenomicRegion reg = GenomicRegion.parse(split[2]);
            Set<Integer> calledBy = new HashSet<>();
            for (int i = 3; i < split.length-2; i++) {
                int called = Integer.parseInt(split[i]);
                if (called == 1) {
                    calledBy.add(i-3);
                } else if (called == 0) {
                    // Not called
                } else {
                    throw new IllegalStateException();
                }
            }
            Tiss tiss = new Tiss(calledBy, Integer.parseInt(split[1]));
            Tsr tsr = reg2tsr.get(reg);
            if (tsr == null) {
                tsr = new Tsr(tiss);
            } else {
                tsr.add(tiss);
            }
            reg2tsr.put(reg, tsr);
        });
        Map<ReferenceSequence, List<Tsr>> outMap = new HashMap<>();
        EI.wrap(tsrMap.keySet()).forEachRemaining(ref -> {
            EI.wrap(tsrMap.get(ref).keySet()).forEachRemaining(reg -> {
                List<Tsr> tsrList = outMap.computeIfAbsent(ref, empty -> new ArrayList<>());
                tsrList.add(tsrMap.get(ref).get(reg));
            });
        });
        return outMap;
    }

    // NOTE: using this will loose some information as the TSR file does not contain information about
    // single Tiss. In this case, use the final TiSS file and the following function "extractTsrsFromFinalFile(...)"
    public static Map<ReferenceSequence, List<Tsr>> extractTsrsFromFinalTsrFile(String path, int skip) throws IOException {
        Map<ReferenceSequence, List<Tsr>> tsrMap = new HashMap<>();
        EI.lines(path).skip(skip).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            List<Tsr> tsrList = tsrMap.computeIfAbsent(Chromosome.obtain(split[0]), empty -> new ArrayList<>());
            GenomicRegion reg = GenomicRegion.parse(split[2]);
            Tsr tsr = new Tsr(reg);
            Set<Integer> calledBy = new HashSet<>();
            for (int i = 3; i < split.length-2; i++) {
                int called = Integer.parseInt(split[i]);
                if (called == 1) {
                    calledBy.add(i-3);
                } else if (called == 0) {
                    // Not called
                } else {
                    throw new IllegalStateException();
                }
            }
            Tiss tiss = new Tiss(calledBy, Integer.parseInt(split[1]));
            tsr.add(tiss);
            tsrList.add(tsr);
        });
        return tsrMap;
    }

    // ========
    // Those two functions are the latest version. use them!
    // ========

    // NOTE: using this will loose some information as the TSR file does not contain information about
    // single Tiss. In this case, use the final TiSS file and the following function "extractTsrsFromFinalTiSSFile(...)"
    public static MemoryIntervalTreeStorage<TsrData> loadTsrsFromFinalTsrFile(String path, int skip) throws IOException {
        MemoryIntervalTreeStorage<TsrData> storage = new MemoryIntervalTreeStorage<>(TsrData.class);
        List<ImmutableReferenceGenomicRegion<TsrData>> entries = new ArrayList<>();
        EI.lines(path).skip(skip).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            ReferenceSequence ref = Chromosome.obtain(split[0]);
            GenomicRegion reg = GenomicRegion.parse(split[2]);
            Set<Integer> calledBy = new HashSet<>();
            for (int i = 3; i < split.length-2; i++) {
                int called = Integer.parseInt(split[i]);
                if (called == 1) {
                    calledBy.add(i-3);
                } else if (called == 0) {
                    // Not called
                } else {
                    throw new IllegalStateException();
                }
            }
            int maxTissPos = Integer.parseInt(split[1]);
            entries.add(new ImmutableReferenceGenomicRegion<>(ref, reg, new TsrData(calledBy.size(), maxTissPos, calledBy)));
        });
        storage.fill(EI.wrap(entries));
        return storage;
    }

    // TsrData will not have the MaxTiSS set!!!
    public static MemoryIntervalTreeStorage<TsrData> loadTsrsFromFinalTissFile(String path, int skip) throws IOException {
        MemoryIntervalTreeStorage<TsrData> storage = new MemoryIntervalTreeStorage<>(TsrData.class);
        Map<ReferenceSequence, Map<GenomicRegion, TsrData>> tsrmap = new HashMap<>();
        EI.lines(path).skip(skip).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            ReferenceSequence ref = Chromosome.obtain(split[0]);
            Map<GenomicRegion, TsrData> regMap = tsrmap.computeIfAbsent(ref, empty -> new HashMap<>());
            GenomicRegion reg = GenomicRegion.parse(split[2]);
            TsrData td = regMap.computeIfAbsent(reg, empty -> new TsrData());
            Set<Integer> calledBy = td.getCalledBy();
            Set<Integer> tissPos = td.getAllTissPos();
            td.setScore(Integer.parseInt(split[split.length-1]));
            for (int i = 3; i < split.length-2; i++) {
                int called = Integer.parseInt(split[i]);
                if (called == 1) {
                    calledBy.add(i-3);
                } else if (called == 0) {
                    // Not called
                } else {
                    throw new IllegalStateException();
                }
            }
            tissPos.add(Integer.parseInt(split[1]));
        });
        List<ImmutableReferenceGenomicRegion<TsrData>> entries = new ArrayList<>();
        EI.wrap(tsrmap.keySet()).forEachRemaining(ref -> {
            EI.wrap(tsrmap.get(ref).keySet()).forEachRemaining(gr -> {
                entries.add(new ImmutableReferenceGenomicRegion<>(ref, gr, tsrmap.get(ref).get(gr)));
            });
        });
        storage.fill(EI.wrap(entries));
        return storage;
    }

    public static void mergeTsrs(Map<ReferenceSequence, List<Tsr>> data) {
        EI.wrap(data.keySet()).forEachRemaining(ref -> {
            List<Tsr> tsrList = data.get(ref);
            tsrList.sort(Tsr::compare);
            if (tsrList.size() < 1) {
                return;
            }
            Tsr lasTsr = tsrList.get(0);
            List<Tsr> toDelete = new ArrayList<>();
            for (int i = 1; i < tsrList.size(); i++) {
                Tsr currentTsr = tsrList.get(i);
                if (currentTsr.getStart() <= lasTsr.getRightmostTissPosition()) {
                    lasTsr.addAll(currentTsr.getMaxTiss());
                    int delta = currentTsr.getEnd() - lasTsr.getEnd();
                    if (delta > 0) {
                        lasTsr.extendBack(delta);
                    }
                    currentTsr.getMaxTiss().clear();
                    toDelete.add(currentTsr);
                } else {
                    lasTsr = currentTsr;
                }
            }
            tsrList.removeAll(toDelete);
        });
    }

    private static Tiss createTiss(int id, int pos) {
        Set<Integer> calledBy = new HashSet<>();
        calledBy.add(id);
        return new Tiss(calledBy, pos);
    }

    public static double log2(double x) {
        return Math.log(x) / Math.log(2) + 1e-10;
    }

    public static double calculateThreshold(List<PeakAndPos> peaks, int avg) {
        if (peaks.size() <= avg) {
            return -9999;
        }
        double xStep = 1.0/peaks.size();
        double yMax = peaks.stream().mapToDouble(PeakAndPos::getValue).max().getAsDouble();
        peaks.sort(Comparator.comparingDouble(PeakAndPos::getValue));
        List<PeakAndPos> peaksNorm = peaks.stream().map(a -> new PeakAndPos(a.getPos(), a.getValue()/yMax)).collect(Collectors.toList());
        peaksNorm.sort(Comparator.comparingDouble(PeakAndPos::getValue));
        double lastMean = lstAvg(peaksNorm, peaksNorm.size() - avg, peaksNorm.size());
        for (int i = peaksNorm.size()-1; i > avg; i--) {
            double currentMean = lstAvg(peaksNorm, i-avg, i);
            if (currentMean - lastMean < xStep) {
                return peaks.get(i-avg/2).getValue();
            }
        }
        return -9999;
    }

    private static double lstAvg(List<PeakAndPos> lst, int from, int to) {
        double sum = 0;
        for (int i = from; i < to; i++) {
            sum += lst.get(i).getValue();
        }
        return sum/(to-from);
    }

    public static int[][] condIndexToCitCondIndex(int[] readConditions, int[] condIndex) {
        int[][] out = new int[readConditions.length][];
        int condAdd = 0;
        for (int i = 0; i < readConditions.length; i++) {
            List<Integer> indexLst = new ArrayList<>();
            int numCond = readConditions[i];
            for (int ci : condIndex) {
                if (ci < numCond + condAdd && ci >= condAdd) {
                    indexLst.add(ci - condAdd);
                }
            }
            condAdd += numCond;
            out[i] = indexLst.stream().mapToInt(Integer::intValue).toArray();
        }
        return out;
    }

    public static float[][] totalsToCitCondIndex(int[] readConditions, float[] totals) {
        float[][] out = new float[readConditions.length][];
        int condAdd = 0;
        for (int i = 0; i < readConditions.length; i++) {
            int numCond = readConditions[i];
            float[] totalTmp = new float[numCond];
            for (int j = 0; j < numCond; j++) {
                totalTmp[j] = totals[condAdd+j];
            }
            condAdd += numCond;
            out[i] = totalTmp;
        }
        return out;
    }

    public static NumericArray extractReadDensities(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                    int[] condIndex, ReferenceSequence ref, int refLength, Strandness strandness) {
        NumericArray readCounts = NumericArray.createMemory(refLength, NumericArray.NumericArrayType.Float);
        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractReadDensitiesFromSingleFile(readCounts, reads.get(citIndex), indices, ref, strandness, new ArrayGenomicRegion(0, refLength), ReadCountMode.Weight);
        }
        return readCounts;
    }

    public static NumericArray extractReadDensities(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                    int[] condIndex, ReferenceSequence ref, int refLength, Strandness strandness, ReadCountMode mode) {
        NumericArray readCounts = NumericArray.createMemory(refLength, NumericArray.NumericArrayType.Float);
        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractReadDensitiesFromSingleFile(readCounts, reads.get(citIndex), indices, ref, strandness, new ArrayGenomicRegion(0, refLength), mode);
        }
        return readCounts;
    }

    public static NumericArray extractReadDensities(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                    int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness) {
        NumericArray readCounts = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractReadDensitiesFromSingleFile(readCounts, reads.get(citIndex), indices, ref, strandness, region, ReadCountMode.Weight);
        }
        return readCounts;
    }

    public static NumericArray extractReadDensitiesNormalized(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                              int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness,
                                                              float[] totals) {
        return extractReadDensitiesNormalized(reads, condIndex, ref, region, strandness, totals, ReadCountMode.Weight);
    }

    public static NumericArray extractReadDensitiesNormalized(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                              int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness,
                                                              float[] totals, ReadCountMode mode) {
        NumericArray readCounts = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        float[][] totalsCitIndices = totalsToCitCondIndex(readsNumConds, totals);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractReadDensitiesFromSingleFileNormalized(readCounts, reads.get(citIndex), indices, ref, strandness, region, totalsCitIndices[citIndex], mode);
        }
        return readCounts;
    }

    public static void extractReadDensitiesFromSingleFile(NumericArray ary, GenomicRegionStorage<AlignedReadsData> reads,
                                                    int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadCountMode mode) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            GenomicRegion reg = r.getRegion();
            double[] counts = r.getData().getTotalCountsForConditions(mode);
            reg.iterator().forEachRemaining(regPart -> {
                for (int i = regPart.getStart(); i < regPart.getEnd(); i++) {
                    int posIndex = i-region.getStart();
                    if (posIndex < 0 || posIndex >= region.getTotalLength()) {
                        continue;
                    }
                    for (int c : condIndex) {
                        ary.setFloat(posIndex, ary.getFloat(posIndex) + (float)counts[c]);
                    }
                }
            });
        });
    }

    public static double[][] extractReadCountsForConditionsFromSingleFile(GenomicRegionStorage<AlignedReadsData> reads,
                                                                          int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadType readType) {
        return extractReadCountsForConditionsFromSingleFile(reads, condIndex, ref, strandness, region, readType, ReadCountMode.Weight);
    }

    public static double[][] extractReadCountsForConditionsFromSingleFile(GenomicRegionStorage<AlignedReadsData> reads,
                                                                             int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadType readType, ReadCountMode readCountMode) {
        if (readType == ReadType.DENSITY) {
            throw new IllegalArgumentException("Only FIVE_PRIME and THREE_PRIME supported. Use extractReadDensitiesForConditionsFromSingleFile() instead.");
        }
        if (strandness == Strandness.Unspecific) {
            throw new IllegalArgumentException("Only Sense and Antisense supported");
        }
        double[][] readCounts = new double[condIndex.length][region.getTotalLength()];
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        final boolean useThreePrimeEnd = strandness == Strandness.Antisense && readType == ReadType.FIVE_PRIME || strandness == Strandness.Sense && readType == ReadType.THREE_PRIME;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            int pos = useThreePrimeEnd ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r);
            double[] counts = r.getData().getTotalCountsForConditions(readCountMode);
            int overlappingRegion = -1;
            for (int or = 0; or < region.getNumParts(); or++) {
                if (region.getPart(or).asRegion().contains(pos)) {
                    overlappingRegion = or;
                    break;
                }
            }
            if (overlappingRegion == -1) {
                return;
            }
            int lengthOfPartsBefore = 0;
            for (int or = 0; or < overlappingRegion; or++) {
                lengthOfPartsBefore += region.getPart(or).asRegion().getTotalLength();
            }
            int posIndex = lengthOfPartsBefore + (pos-region.getPart(overlappingRegion).getStart());
//            if (posIndex < 0 || posIndex >= region.getTotalLength()) {
//                return;
//            }
            for (int c = 0; c < condIndex.length; c++) {
                readCounts[c][posIndex] += counts[condIndex[c]];
            }
        });

        return readCounts;
    }

    public static double[][] extractReadDensitiesForConditionsFromSingleFile(GenomicRegionStorage<AlignedReadsData> reads,
                                                          int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadCountMode readCountMode) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : strandness.equals(Strandness.Unspecific) ? ref.toStrandIndependent() : ref;
        double[][] readDensities = new double[condIndex.length][region.getTotalLength()];
        reads.ei(refTmp, region).forEachRemaining(r -> {
            GenomicRegion reg = r.getRegion();
            double[] counts = r.getData().getTotalCountsForConditions(readCountMode);
            reg.iterator().forEachRemaining(regPart -> {
                for (int i = regPart.getStart(); i < regPart.getEnd(); i++) {
                    int overlappingRegion = -1;
                    for (int or = 0; or < region.getNumParts(); or++) {
                        if (region.getPart(or).asRegion().contains(i)) {
                            overlappingRegion = or;
                            break;
                        }
                    }
                    if (overlappingRegion == -1) {
                        continue;
                    }
                    int lengthOfPartsBefore = 0;
                    for (int or = 0; or < overlappingRegion; or++) {
                        lengthOfPartsBefore += region.getPart(or).asRegion().getTotalLength();
                    }
                    int posIndex = lengthOfPartsBefore + (i-region.getPart(overlappingRegion).getStart());
                    for (int c = 0; c < condIndex.length; c++) {
                        readDensities[c][posIndex] += counts[condIndex[c]];
                    }
                }
            });
        });
        return readDensities;
    }

    public static void extractReadDensitiesFromSingleFileNormalized(NumericArray ary, GenomicRegionStorage<AlignedReadsData> reads,
                                                                    int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region,
                                                                    float[] totals) {
        extractReadDensitiesFromSingleFileNormalized(ary, reads, condIndex, ref, strandness, region, totals, ReadCountMode.Weight);
    }

    public static void extractReadDensitiesFromSingleFileNormalized(NumericArray ary, GenomicRegionStorage<AlignedReadsData> reads,
                                                                    int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region,
                                                                    float[] totals, ReadCountMode mode) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            GenomicRegion reg = r.getRegion();
            double[] counts = r.getData().getTotalCountsForConditions(mode);
            reg.iterator().forEachRemaining(regPart -> {
                for (int i = regPart.getStart(); i < regPart.getEnd(); i++) {
                    int posIndex = i-region.getStart();
                    if (posIndex < 0 || posIndex >= region.getTotalLength()) {
                        continue;
                    }
                    for (int c : condIndex) {
                        ary.setFloat(posIndex, ary.getFloat(posIndex) + (((float)counts[c])/totals[c])*1000000);
                    }
                }
            });
        });
    }

    public static double getMappability(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                        int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness) {
        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        double readsAll = 0;
        double readsWeighted = 0;
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            Mappability m = getMappabilityFromSingleFile(reads.get(citIndex), indices, ref, strandness, region);
            readsAll += m.readSumAll;
            readsWeighted += m.readSumWeighted;
        }
        return readsWeighted / readsAll;
    }

    public static Mappability getMappabilityFromSingleFile(GenomicRegionStorage<AlignedReadsData> reads,
                                                      int[] condIndex, ReferenceSequence ref, Strandness strandness, GenomicRegion region) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        double readAll = 0;
        double readWeight = 0;
        for (ReferenceGenomicRegion<AlignedReadsData> r : reads.ei(refTmp, region).loop()) {
            GenomicRegion reg = r.getRegion();
            double[] countsWeighted = r.getData().getTotalCountsForConditions(ReadCountMode.Weight);
            double[] countsAll = r.getData().getTotalCountsForConditions(ReadCountMode.All);
            for (GenomicRegionPart regPart : reg.iterator().loop()) {
                for (int i = regPart.getStart(); i < regPart.getEnd(); i++) {
                    int posIndex = i-region.getStart();
                    if (posIndex < 0 || posIndex >= region.getTotalLength()) {
                        continue;
                    }
                    for (int c : condIndex) {
                        readWeight += countsWeighted[c];
                        readAll += countsAll[c];
                    }
                }
            }
        }
        return new Mappability(readWeight, readAll);
    }

    //TODO: put somewhere else!
    public static class Mappability {
        public double readSumWeighted;
        public double readSumAll;

        public Mappability(double readSumWeighted, double readSumAll) {
            this.readSumWeighted = readSumWeighted;
            this.readSumAll = readSumAll;
        }
    }

    public static NumericArray extractFivePrimeCounts(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                      int[] condIndex, ReferenceSequence ref, int refLength, Strandness strandness) {

        NumericArray readCounts = NumericArray.createMemory(refLength, NumericArray.NumericArrayType.Float);
        for (int cond : condIndex) {
            int lstIndex = 0;
            while (cond >= reads.get(lstIndex).getMetaDataConditions().length) {
                cond -= reads.get(lstIndex).getMetaDataConditions().length;
                lstIndex++;
            }
            extractCountsFromSingleFile(readCounts, reads.get(lstIndex), cond, ref, strandness, new ArrayGenomicRegion(0, refLength), ReadType.FIVE_PRIME);
        }
        return readCounts;
    }

    public static NumericArray extractCounts(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                      int[] condIndex, ReferenceSequence ref, int refLength, Strandness strandness, ReadType readType) {

        NumericArray readCounts = NumericArray.createMemory(refLength, NumericArray.NumericArrayType.Float);
        for (int cond : condIndex) {
            int lstIndex = 0;
            while (cond >= reads.get(lstIndex).getMetaDataConditions().length) {
                cond -= reads.get(lstIndex).getMetaDataConditions().length;
                lstIndex++;
            }
            extractCountsFromSingleFile(readCounts, reads.get(lstIndex), cond, ref, strandness, new ArrayGenomicRegion(0, refLength), readType);
        }
        return readCounts;
    }

    public static NumericArray extractFivePrimeCounts(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                      int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness) {

        NumericArray readCounts = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
        for (int cond : condIndex) {
            int lstIndex = 0;
            while (cond >= reads.get(lstIndex).getMetaDataConditions().length) {
                cond -= reads.get(lstIndex).getMetaDataConditions().length;
                lstIndex++;
            }
            extractCountsFromSingleFile(readCounts, reads.get(lstIndex), cond, ref, strandness, region, ReadType.FIVE_PRIME);
        }
        return readCounts;
    }

    public static double[] extractShortTotalRatios(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                          int[] condIndex, int maxShortLength, ReferenceSequence ref, GenomicRegion region, Strandness strandness) {
        double[] shortTotalRatios = new double[region.getTotalLength()];

        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractShortTotalRatiosFromSingleFile(shortTotalRatios, reads.get(citIndex), indices, maxShortLength, ref, strandness, region, ReadCountMode.Weight);
        }

        return shortTotalRatios;
    }

    public static void extractShortTotalRatiosFromSingleFile(double[] shortTotalRatios, GenomicRegionStorage<AlignedReadsData> reads, int[] conds, int maxShortLength,
                                                             ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadCountMode mode) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        double[][] shortAndTotalReads = new double[2][shortTotalRatios.length];
        reads.ei(refTmp, region).filter(r -> region.contains(strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r))).forEachRemaining(r -> {
            GenomicRegion reg = r.getRegion();
            double[] counts = r.getData().getTotalCountsForConditions(mode);
            int readLength = reg.getEnd() - reg.getStart();
            int pos = (strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r)) - region.getStart();
            double count = 0.0;
            for (int c : conds) {
                count += counts[c];
            }
            if (readLength <= maxShortLength) {
                shortAndTotalReads[0][pos] += count;
            }
            shortAndTotalReads[1][pos] += count;
        });
        for (int i = 0; i < shortTotalRatios.length; i++) {
            shortTotalRatios[i] = shortAndTotalReads[1][i] == 0 ? 0 : shortAndTotalReads[0][i]/shortAndTotalReads[1][i];
        }
    }

    public static Map<Integer, Double> extractReadLengths(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                           int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness) {
        Map<Integer, Double> readLengths = new HashMap<>();

        int[] readsNumConds = reads.stream().mapToInt(r -> r.getMetaDataConditions().length).toArray();
        int[][] condCitIndices = condIndexToCitCondIndex(readsNumConds, condIndex);
        for (int citIndex = 0; citIndex < condCitIndices.length; citIndex++) {
            int[] indices = condCitIndices[citIndex];
            extractReadLengthsFromSingleFile(readLengths, reads.get(citIndex), indices, ref, strandness, region, ReadCountMode.Weight);
        }

        return readLengths;
    }

    public static void extractReadLengthsFromSingleFile(Map<Integer, Double> readLengths, GenomicRegionStorage<AlignedReadsData> reads, int[] conds,
                                                        ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadCountMode mode) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        reads.ei(refTmp, region).filter(r -> region.contains(strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r))).forEachRemaining(r -> {
            GenomicRegion reg = r.getRegion();
            double[] counts = r.getData().getTotalCountsForConditions(mode);
            int readLength = reg.getEnd() - reg.getStart();
            if (!readLengths.containsKey(readLength)) {
                readLengths.put(readLength, 0.0);
            }
            double count = 0.0;
            for (int c : conds) {
                count += counts[c];
            }
            readLengths.put(readLength, readLengths.get(readLength) + count);
        });
    }

    /**
     * normalized based on read-totals. REGION LENGTH IS NOT CONSIDERED, SO NO RPKM!!!!
     * @param reads
     * @param condIndex
     * @param ref
     * @param region
     * @param strandness
     * @param totals needs to be the same length as conditions in reads
     * @return
     */
    public static NumericArray extractFivePrimeCountsNormalized(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                                int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness,
                                                                float[] totals) {

        NumericArray readCounts = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
        for (int cond : condIndex) {
            int lstIndex = 0;
            while (cond >= reads.get(lstIndex).getMetaDataConditions().length) {
                cond -= reads.get(lstIndex).getMetaDataConditions().length;
                lstIndex++;
            }
            NumericArray readCountsTmp = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
            extractCountsFromSingleFile(readCountsTmp, reads.get(lstIndex), cond, ref, strandness, region, ReadType.FIVE_PRIME);
            for (int i = 0; i < readCountsTmp.length(); i++) {
                readCounts.setFloat(i, readCounts.getFloat(i) + ((readCountsTmp.getFloat(i)/totals[cond]) * 1.0E6f));
            }
        }
        return readCounts;
    }

    public static void extractFivePrimeCountsFromSingleFile(NumericArray ary, GenomicRegionStorage<AlignedReadsData> reads, int cond, ReferenceSequence ref, Strandness strandness) {
        extractCountsFromSingleFile(ary, reads, cond, ref, strandness, new ArrayGenomicRegion(0, ary.length()), ReadType.FIVE_PRIME);
    }

    public static void extractCountsFromSingleFile(NumericArray ary, GenomicRegionStorage<AlignedReadsData> reads, int cond, ReferenceSequence ref, Strandness strandness, GenomicRegion region, ReadType readType) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        final boolean switchFiveAndThreePrimeEnd = strandness == Strandness.Antisense && readType == ReadType.FIVE_PRIME || strandness == Strandness.Sense && readType == ReadType.THREE_PRIME;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            int pos0 = switchFiveAndThreePrimeEnd ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r);
            int pos1 = switchFiveAndThreePrimeEnd ? GenomicRegionPosition.ThreePrime.position(r, 1) : GenomicRegionPosition.FivePrime.position(r, 1);
            NumericArray c0 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            NumericArray c1 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            for (int k = 0; k < r.getData().getDistinctSequences(); k++) {
                if (hasEndMismatch(r.getData(), k, r.getRegion().getTotalLength(), strandness)) {
                    c1 = r.getData().addCountsForDistinct(k, c1, ReadCountMode.Weight);
                } else {
                    c0 = r.getData().addCountsForDistinct(k, c0, ReadCountMode.Weight);
                }
            }
            int posIndex0 = pos0 - region.getStart();
            if (posIndex0 >= 0 && posIndex0 < ary.length()) {
                if (c0.length() > 0) {
                    ary.setFloat(posIndex0, ary.getFloat(posIndex0) + c0.getFloat(cond));
                }
            }
            int posIndex1 = pos1 - region.getStart();
            if (posIndex1 >= 0 && posIndex1 < ary.length()) {
                if (c1.length() > 0) {
                    ary.setFloat(posIndex1, ary.getFloat(posIndex1) + c1.getFloat(cond));
                }
            }
        });
    }

    public static SparseNumericArray<Double> extractFivePrimeCountsNormalizedSparse(List<GenomicRegionStorage<AlignedReadsData>> reads,
                                                                int[] condIndex, ReferenceSequence ref, GenomicRegion region, Strandness strandness,
                                                                float[] totals) {
        SparseNumericArray<Double> readCounts = new SparseNumericArray<>(region.getTotalLength(), 0.0d);
        for (int cond : condIndex) {
            int lstIndex = 0;
            while (cond >= reads.get(lstIndex).getMetaDataConditions().length) {
                cond -= reads.get(lstIndex).getMetaDataConditions().length;
                lstIndex++;
            }
            SparseNumericArray<Double> readCountsTmp = new SparseNumericArray<>(region.getTotalLength(), 0.0d);
            extractFivePrimeCountsFromSingleFileSparse(readCountsTmp, reads.get(lstIndex), cond, ref, strandness, region);
            for (int i : readCountsTmp.getNonZeroIndices()) {
                readCounts.set(i, readCounts.get(i) + ((readCountsTmp.get(i)/totals[cond]) * 1.0E6d));
            }
        }
        return readCounts;
    }

    public static void extractFivePrimeCountsFromSingleFileSparse(SparseNumericArray<Double> ary, GenomicRegionStorage<AlignedReadsData> reads, int cond, ReferenceSequence ref, Strandness strandness) {
        extractFivePrimeCountsFromSingleFileSparse(ary, reads, cond, ref, strandness, new ArrayGenomicRegion(0, ary.length()));
    }

    public static void extractFivePrimeCountsFromSingleFileSparse(SparseNumericArray<Double> ary, GenomicRegionStorage<AlignedReadsData> reads, int cond, ReferenceSequence ref, Strandness strandness, GenomicRegion region) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            int pos0 = strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r);
            int pos1 = strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r, 1) : GenomicRegionPosition.FivePrime.position(r, 1);
            NumericArray c0 = NumericArray.createMemory(0, NumericArray.NumericArrayType.Double);
            NumericArray c1 = NumericArray.createMemory(0, NumericArray.NumericArrayType.Double);
            for (int k = 0; k < r.getData().getDistinctSequences(); k++) {
                if (hasEndMismatch(r.getData(), k, r.getRegion().getTotalLength(), strandness)) {
                    c1 = r.getData().addCountsForDistinct(k, c1, ReadCountMode.Weight);
                } else {
                    c0 = r.getData().addCountsForDistinct(k, c0, ReadCountMode.Weight);
                }
            }
            int posIndex0 = pos0 - region.getStart();
            if (posIndex0 >= 0 && posIndex0 < ary.length()) {
                if (c0.length() > 0) {
                    ary.set(posIndex0, ary.get(posIndex0) + c0.getDouble(cond));
                }
            }
            int posIndex1 = pos1 - region.getStart();
            if (posIndex1 >= 0 && posIndex1 < ary.length()) {
                if (c1.length() > 0) {
                    ary.set(posIndex1, ary.get(posIndex1) + c1.getDouble(cond));
                }
            }
        });
    }

    public static void extractFivePrimeCountsFromSingleFileSparseFull(SparseNumericArray<Float>[] ary, GenomicRegionStorage<AlignedReadsData> reads, ReferenceSequence ref, Strandness strandness) {
        extractFivePrimeCountsFromSingleFileSparseFull(ary, reads, ref, strandness, new ArrayGenomicRegion(0, ary[0].length()));
    }

    public static void extractFivePrimeCountsFromSingleFileSparseFull(SparseNumericArray<Float>[] ary, GenomicRegionStorage<AlignedReadsData> reads, ReferenceSequence ref, Strandness strandness, GenomicRegion region) {
        ReferenceSequence refTmp = strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref;
        reads.ei(refTmp, region).forEachRemaining(r -> {
            int pos0 = strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r);
            int pos1 = strandness.equals(Strandness.Antisense) ? GenomicRegionPosition.ThreePrime.position(r, 1) : GenomicRegionPosition.FivePrime.position(r, 1);
            NumericArray c0 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            NumericArray c1 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            for (int k = 0; k < r.getData().getDistinctSequences(); k++) {
                if (hasEndMismatch(r.getData(), k, r.getRegion().getTotalLength(), strandness)) {
                    c1 = r.getData().addCountsForDistinct(k, c1, ReadCountMode.Weight);
                } else {
                    c0 = r.getData().addCountsForDistinct(k, c0, ReadCountMode.Weight);
                }
            }
            int posIndex0 = pos0 - region.getStart();
            if (posIndex0 >= 0 && posIndex0 < ary[0].length()) {
                if (c0.length() > 0) {
                    for (int i = 0; i < c0.length(); i++) {
                        ary[i].set(posIndex0, ary[i].get(posIndex0) + c0.getFloat(i));
                    }
                }
            }
            int posIndex1 = pos1 - region.getStart();
            if (posIndex1 >= 0 && posIndex1 < ary[0].length()) {
                if (c1.length() > 0) {
                    for (int i = 0; i < c1.length(); i++) {
                        ary[i].set(posIndex1, ary[i].get(posIndex1) + c1.getFloat(i));
                    }
                }
            }
            if(c0.length() > 0 && c0.length() != ary.length || c1.length() > 0 && c1.length() != ary.length) {
                System.err.println("Wrong lengths");
            }
        });
    }

    private static boolean hasEndMismatch(AlignedReadsData ard, int distict, int readLength, Strandness strandness) {
        if (strandness.equals(Strandness.Antisense)) {
            return hasTailingMismatch(ard, distict, readLength);
        } else {
            return RiboUtils.hasLeadingMismatch(ard, distict);
        }
    }

    private static boolean hasTailingMismatch(AlignedReadsData ard, int distinct, int readLength) {
        for (int i=0; i<ard.getVariationCount(distinct); i++) {
            if (ard.isMismatch(distinct, i) && ard.getMismatchPos(distinct, i)==readLength-1)
                return true;
        }
        return false;
    }

    public static GenomicRegion createRegion(int pos, int windowSize, boolean toTheRight) {
        int regStart = (toTheRight ? pos : pos - windowSize);
        int regStop = (toTheRight ? pos + windowSize: pos);
        return new ArrayGenomicRegion(regStart, regStop);
    }

    public static List<MutableTriple<Integer, Double, Double>> cleanUpMultiValueData(List<MutableTriple<Integer, Double, Double>> lstToClean, int multiThreshold) {
        Map<Long, Integer> countMap = new HashMap<>();
        for (MutableTriple<Integer, Double, Double> itm : lstToClean) {
            long bits = Double.doubleToLongBits(itm.Item2);
            if (countMap.containsKey(bits)) {
                countMap.put(bits, countMap.get(bits) + 1);
            } else {
                countMap.put(bits, 1);
            }
        }
        return EI.wrap(lstToClean).filter(i -> countMap.get(Double.doubleToLongBits(i.Item2)) < multiThreshold).list();
    }

    public static List<MutablePair<Integer, Double>> cleanUpMultiValueDataPair(List<MutablePair<Integer, Double>> lstToClean, int multiThreshold) {
        Map<Long, Integer> countMap = new HashMap<>();
        for (MutablePair<Integer, Double> itm : lstToClean) {
            long bits = Double.doubleToLongBits(itm.Item2);
            if (countMap.containsKey(bits)) {
                countMap.put(bits, countMap.get(bits) + 1);
            } else {
                countMap.put(bits, 1);
            }
        }
        return EI.wrap(lstToClean).filter(i -> countMap.get(Double.doubleToLongBits(i.Item2)) < multiThreshold).list();
    }

    public static List<MutableTriple<Integer, Double, Double>> cleanUpMultiValueDataTriple(List<MutableTriple<Integer, Double, Double>> lstToClean, int multiThreshold) {
        Map<Long, Integer> countMap = new HashMap<>();
        for (MutableTriple<Integer, Double, Double> itm : lstToClean) {
            long bits = Double.doubleToLongBits(itm.Item2);
            if (countMap.containsKey(bits)) {
                countMap.put(bits, countMap.get(bits) + 1);
            } else {
                countMap.put(bits, 1);
            }
        }
        return EI.wrap(lstToClean).filter(i -> countMap.get(Double.doubleToLongBits(i.Item2)) < multiThreshold).list();
    }

    public static void addNeighbourWeights(List<MutableTriple<Integer, Double, Double>> sortedLst) {
        for (int i = 0; i < sortedLst.size(); i++) {
            MutableTriple<Integer, Double, Double> current = sortedLst.get(i);

            int j = i-1;
            while (j >= 0 && current.Item1 - sortedLst.get(j).Item1 < 100) {
                int dist = current.Item1 - sortedLst.get(j).Item1;
                double weight = (double)dist/100.d;
                current.Item3 += sortedLst.get(j).Item2 * weight;
                j--;
            }

            j = i+1;
            while (j < sortedLst.size() && sortedLst.get(j).Item1 - current.Item1 < 100) {
                int dist = sortedLst.get(j).Item1 - current.Item1;
                double weight = (double)dist/100.d;
                current.Item3 += sortedLst.get(j).Item2 * weight;
                j++;
            }
        }
    }

    public static void mergeThresholdData(String[] inFiles, String outFile, int window, boolean isPVal) throws IOException {
//        mergeCustomData(inFiles, outFile, window, isPVal, 1, 0, -1, 1, 2, true, 100);
        Map<ReferenceSequence, List<MutableTriple<Integer, Double, Double>>> map = new HashMap<>();
        for (String inFile : inFiles) {
            EI.lines(inFile).skip(1).forEachRemaining(l -> {
                String[] split = StringUtils.split(l, "\t");
                ReferenceSequence ref = Chromosome.obtain(split[0]);
                int pos = Integer.parseInt(split[1]);
                double val;
                if (split.length > 3) {
                    double v1 = Double.parseDouble(split[2]);
                    double v2 = Double.parseDouble(split[2+1]);
                    val = v1 > v2 ? v1 : v2;
                } else if (isPVal) {
                    val = 1 - Double.parseDouble(split[2]);
                } else {
                    val = Double.parseDouble(split[2]);
                }
                map.computeIfAbsent(ref, absent -> new ArrayList<>()).add(new MutableTriple<>(pos, val, val));
            });
        }
        LineWriter writer = new LineOrientedFile(outFile).write();
        writer.writeLine("Ref\tTiSS\tValue");
        EI.wrap(map.keySet()).forEachRemaining(ref -> {
            List<MutableTriple<Integer, Double, Double>> lst = map.get(ref);
            lst.sort(Comparator.comparingInt(a -> a.Item1));
            addNeighbourWeights(lst);
            lst = cleanUpMultiValueData(lst, 100);
            int lastPos = lst.get(0).Item1;
            MutableTriple<Integer, Double, Double> lastHighest = lst.get(0);
            for (int i = 1; i < lst.size(); i++) {
                MutableTriple<Integer, Double, Double> current = lst.get(i);
                if (current.Item1-lastPos > window) {
                    try {
                        writer.writeLine(ref + "\t" + lastHighest.Item1 + "\t" + lastHighest.Item3);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    lastHighest = current;
                } else {
                    lastHighest = current.Item3 > lastHighest.Item3 ? current : lastHighest;
                }
                lastPos = current.Item1;
            }
            try {
                writer.writeLine(ref + "\t" + lastHighest.Item1 + "\t" + lastHighest.Item3);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        writer.close();
    }

    public static void normalize(Map<ReferenceSequence, List<MutablePair<Integer, Double>>> map) {
        double maxVal = -9999, minVal = 9999999;
        for (List<MutablePair<Integer, Double>> lst : map.values()) {
            for (MutablePair<Integer, Double> p : lst) {
                maxVal = p.Item2 > maxVal ? p.Item2 : maxVal;
                minVal = p.Item2 < minVal ? p.Item2 : minVal;
            }
        }
        maxVal -= minVal;
        for (List<MutablePair<Integer, Double>> lst : map.values()) {
            for (MutablePair<Integer, Double> p : lst) {
                p.Item2 = (p.Item2-minVal) / maxVal;
            }
        }
    }

    public static String getGeneNameFromGeneId(String geneId, Genomic genomic) {
        return genomic.getGeneTable("symbol").apply(geneId);
    }

    public static ImmutableReferenceGenomicRegion<NameAnnotation> getRgrFromGeneName(String geneName, Genomic genomic) {
        ReferenceGenomicRegion<?> rgr = genomic.getNameIndex().get(geneName);
        if (rgr == null) {
            return null;
        }
        return new ImmutableReferenceGenomicRegion<>(rgr.getReference(), rgr.getRegion(), new NameAnnotation(rgr.getData().toString()));
    }

    public static int[] getOrder(double[] ary, Comparator<Double> doubleComparator) {
        ArrayList<Double> cloned = EI.wrap(ary).list();
        ArrayList<Double> clonedToSort = EI.wrap(ary).list();
        clonedToSort.sort(doubleComparator);
        int skip = 0;
        Set<Integer> indices = EI.seq(0, ary.length).set();
        int[] out = new int[ary.length];
        for (int i = 0; i < out.length; i++) {
            int index = indexOf(cloned, clonedToSort.get(i), skip);
            if (indices.contains(index)) {
                indices.remove(index);
                out[i] = index;
                skip = 0;
            } else {
                skip++;
                i--;
            }
        }
        return out;
    }

    private static int indexOf(List<Double> lst, Double d, int skip) {
        int skipped = 0;
        for (int i = 0; i < lst.size(); i++) {
            if (lst.get(i).equals(d)) {
                if (skipped == skip) {
                    return i;
                }
                skipped++;
            }
        }
        return -1;
    }

    public static double mean(List<Integer> lst) {
        double total = lst.stream().mapToDouble(e -> (double) e).sum();
        return total/lst.size();
    }

    public static double sampleVar(List<Integer> lst) {
        double mean = mean(lst);
        return lst.stream().mapToDouble(e -> Math.pow(e-mean, 2.d)).sum() / (lst.size()-1);
    }

    public static double sampleSd(List<Integer> lst) {
        return Math.sqrt(sampleVar(lst));
    }

    public static double getReadCount(CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> cit, ReferenceSequence ref, GenomicRegion reg, int[] conditions, Strandness strandness) {
        AtomicReference<Double> readCount = new AtomicReference<>(0.0);
        ReferenceSequence refTmp = strandness == Strandness.Antisense ? ref.toOppositeStrand() : ref;
        cit.ei(refTmp, reg).forEachRemaining(r -> {
            double[] counts = r.getData().getTotalCountsForConditions(ReadCountMode.Weight);
            for (int i : conditions) {
                readCount.updateAndGet(v -> v + counts[i]);
            }
        });
        return readCount.get();
    }

    public static double getMaxReadCountFivePrime(GenomicRegionStorage<AlignedReadsData> reads, int cond, ReferenceSequence ref, Strandness strandness, GenomicRegion region) {
        NumericArray ary = NumericArray.createMemory(region.getTotalLength(), NumericArray.NumericArrayType.Float);
        extractCountsFromSingleFile(ary, reads, cond, ref, strandness, region, ReadType.FIVE_PRIME);
        return ArrayUtils.max(ary.toDoubleArray());
    }

    public static double getReadCountFivePrime(CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> cit, ReferenceSequence ref, GenomicRegion reg, int[] conditions, Strandness strandness, ReadCountMode readCountMode) {
        AtomicReference<Double> readCount = new AtomicReference<>(0.0);
        ReferenceSequence refTmp = strandness == Strandness.Antisense ? ref.toOppositeStrand() : ref;
        cit.ei(refTmp, reg).filter(f -> {
            if (strandness == Strandness.Antisense) {
                return reg.contains(GenomicRegionPosition.ThreePrime.position(f));
            } else {
                return reg.contains(GenomicRegionPosition.FivePrime.position(f));
            }
        }).forEachRemaining(r -> {
            double[] counts = r.getData().getTotalCountsForConditions(readCountMode);
            for (int i : conditions) {
                readCount.updateAndGet(v -> v + counts[i]);
            }
        });
        return readCount.get();
    }

    public static double getCountFromRmq(DiskGenomicNumericProvider rmq, ReferenceSequence ref, GenomicRegion reg, final int[] conditions) {
        double[] counts = new double[rmq.getNumDataRows()];
        GenomicNumericProvider.PositionNumericIterator it = rmq.iterateValues(ref, reg);
        while(it.hasNext()) {
            it.nextInt();
            ArrayUtils.add(counts, it.getValues(null));
        }
        double readCount = 0;
        for (int c : conditions) {
            readCount += counts[c];
        }
        return readCount;
    }

    public static double[] getReadCountForConditions(CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> cit, ReferenceSequence ref, GenomicRegion reg, int[] conditions, Strandness strandness, ReadCountMode readCountMode) {
        double[] out= new double[conditions.length];
        ReferenceSequence refTmp = strandness == Strandness.Antisense ? ref.toOppositeStrand() : ref;
        cit.ei(refTmp, reg).forEachRemaining(r -> {
            double[] counts = r.getData().getTotalCountsForConditions(readCountMode);
            for (int i = 0; i < conditions.length; i++) {
                out[i] += counts[conditions[i]];
            }
        });
        return out;
    }

    public static int getDistance(ReferenceGenomicRegion rgr1, ReferenceGenomicRegion rgr2) {
        if (!rgr1.getReference().equals(rgr2.getReference())) {
            return -1;
        }

        if (rgr1.getRegion().intersects(rgr2.getRegion())) {
            return 0;
        }
        int s1 = rgr1.getRegion().getStart();
        int s2 = rgr2.getRegion().getStart();
        int e1 = rgr1.getRegion().getEnd();
        int e2 = rgr2.getRegion().getEnd();
        if (s1 < s2) {
            return s2-e1;
        } else {
            return s1-e2;
        }
    }

    public static double[] getReadCountFivePrimeForConditions(CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> cit, ReferenceSequence ref, GenomicRegion reg, int[] conditions, Strandness strandness, ReadCountMode readCountMode) {
        double[] out= new double[conditions.length];
        ReferenceSequence refTmp = strandness == Strandness.Antisense ? ref.toOppositeStrand() : ref;
        cit.ei(refTmp, reg).filter(f -> {
            if (strandness == Strandness.Antisense) {
                return reg.contains(GenomicRegionPosition.ThreePrime.position(f));
            } else {
                return reg.contains(GenomicRegionPosition.FivePrime.position(f));
            }
        }).forEachRemaining(r -> {
            double[] counts = r.getData().getTotalCountsForConditions(readCountMode);
            for (int i = 0; i < conditions.length; i++) {
                out[i] += counts[conditions[i]];
            }
        });
        return out;
    }

    public static AlignedReadsData alignedReadsDataToOppositeStrand(AlignedReadsData data) {
        AlignedReadsDataFactory factory = new AlignedReadsDataFactory(data.getNumConditions(), data.hasNonzeroInformation());

        factory.start(); // Start factory
        for (int d = 0; d < data.getDistinctSequences(); d++) {
            factory.newDistinctSequence(); // add distinct sequence
            if (data.hasGeometry()) { // check if paired-end
                factory.setGeometry(data.getGeometryBeforeOverlap(d), data.getGeometryOverlap(d), data.getGeometryAfterOverlap(d));
            }
            factory.setMultiplicity(data.getMultiplicity(d));
            if (data.hasWeights()) {
                factory.setWeight(d, data.getWeight(d));
            }
            if (data.hasNonzeroInformation()) { // check if sparse-data and set counts accordingly
                for (int i : data.getNonzeroCountIndicesForDistinct(d)) {
                    factory.setCount(i, data.getCount(d, i));
                }
            }
            else {
                for (int c = 0; c < data.getNumConditions(); c++) {
                    int count = data.getCount(d, c);
                    if (count > 0) {
                        factory.setCount(d, c, count);
                    }
                }
            }
            for (AlignedReadsVariation var : data.getVariations(d)) { // convert variations to other strand.
                // Careful: Positions of variations are in 5prime to 3prime orientation, which needs to be swapped.
                if (var.isSoftclip()) {
                    AlignedReadsSoftclip softclip = (AlignedReadsSoftclip) var;
                    factory.addVariation(new AlignedReadsSoftclip(softclip.getPosition()!=0,
                                                                            SequenceUtils.getDnaReverseComplement(softclip.getReadSequence()),
                                                                            softclip.isFromSecondRead()));
                } else if (var.isMismatch()) {
                    AlignedReadsMismatch mismatch = (AlignedReadsMismatch) var;
                    factory.addVariation(new AlignedReadsMismatch(data.getMappedLength(d) - mismatch.getPosition() - 1,
                                                                            SequenceUtils.getDnaReverseComplement(mismatch.getReferenceSequence()),
                                                                            SequenceUtils.getDnaReverseComplement(mismatch.getReadSequence()), mismatch.isFromSecondRead()));
                } else if (var.isDeletion()) {
                    AlignedReadsDeletion deletion = (AlignedReadsDeletion) var;
                    factory.addVariation(new AlignedReadsDeletion(data.getMappedLength(d) - deletion.getPosition() - 1,
                                                                            SequenceUtils.getDnaReverseComplement(deletion.getReferenceSequence()),
                                                                            deletion.isFromSecondRead()));
                } else if (var.isInsertion()) {
                    AlignedReadsInsertion insertion = (AlignedReadsInsertion) var;
                    factory.addVariation(new AlignedReadsInsertion(data.getMappedLength(d) - insertion.getPosition() - 1,
                                                                            insertion.getReadSequence(),
                                                                            insertion.isFromSecondRead()));
                }
            }
            if (data.hasId()) {
                factory.setId(d, data.getId(d));
            }
        }
        factory.makeDistinct();

        return factory.create();
    }

    public static int[][] conditionMergeStringToIntArray(String toMerge) {
        return EI.wrap(StringUtils.split(toMerge, "/")).map(m -> EI.wrap(StringUtils.split(m, ",")).mapToInt(Integer::parseInt).toIntArray()).toArray(new int[0][]);
    }

    public static <T> ImmutableReferenceGenomicRegion<T> getClosestReferenceGenomicRegionDownstream(ReferenceGenomicRegion<?> rgr, GenomicRegionStorage<T> storage, int maxDist) {
        return getClosestReferenceGenomicRegionDownstream(rgr.getReference(), rgr.getRegion(), storage, maxDist);
    }

    public static <T> ImmutableReferenceGenomicRegion<T> getClosestReferenceGenomicRegionDownstream(ReferenceSequence ref, GenomicRegion reg, GenomicRegionStorage<T> storage, int maxDist) {
        ImmutableReferenceGenomicRegion<T> closest = null;

        if (ref.isPlus()) {
            ExtendedIterator<ImmutableReferenceGenomicRegion<T>> ei = storage.ei(ref, new ArrayGenomicRegion(reg.getEnd(), reg.getEnd() + maxDist)).filter(f -> f.getRegion().getStart() > reg.getEnd());
            if (!ei.hasNext()) {
                return closest;
            }
            closest = ei.next();
            for (ImmutableReferenceGenomicRegion<T> current : ei.loop()) {
                if (current.getRegion().getStart() - reg.getEnd() < closest.getRegion().getStart() - reg.getEnd()) {
                    closest = current;
                }
            }
        } else if (ref.isMinus()) {
            ExtendedIterator<ImmutableReferenceGenomicRegion<T>> ei = storage.ei(ref, new ArrayGenomicRegion(reg.getStart() - maxDist, reg.getStart())).filter(f -> f.getRegion().getEnd() < reg.getStart());
            if (!ei.hasNext()) {
                return closest;
            }
            closest = ei.next();
            for (ImmutableReferenceGenomicRegion<T> current : ei.loop()) {
                if (reg.getStart() - current.getRegion().getEnd() < reg.getStart() - closest.getRegion().getEnd()) {
                    closest = current;
                }
            }
        } else {
            throw new UnsupportedOperationException("Only plus or minus references allowed for now.");
        }

        return closest;
    }

    public static <T> List<ImmutableReferenceGenomicRegion<T>> getClosestDownstreamRgrInStorage(ReferenceGenomicRegion<?> rgr, GenomicRegionStorage<T> storage, int maxDist) {
        List<ImmutableReferenceGenomicRegion<T>> closest = new ArrayList<>();

        ArrayList<ImmutableReferenceGenomicRegion<T>> sorted = storage.ei(rgr.getDownstream(maxDist)).
                filter(f -> f.getReference().isPlus() ? f.getRegion().getStart() > rgr.getRegion().getStart() : f.getRegion().getStop() < rgr.getRegion().getStop()).
                sort(Comparator.comparingInt(s -> s.getReference().isPlus() ? s.getRegion().getStart() : -s.getRegion().getStop())).list();

        if (sorted.size() == 0) {
            return closest;
        }

        for (ImmutableReferenceGenomicRegion<T> e : sorted) {
            if (closest.size() == 0) {
                closest.add(e);
                continue;
            }
            if (e.getReference().isPlus() ? e.getRegion().getStart() > closest.get(0).getRegion().getStart() : e.getRegion().getStop() < closest.get(0).getRegion().getStop()) {
                break;
            }
            closest.add(e);
        }

        return closest;
    }

    public static class TtoCData {
        public double totalT;
        public double totalTtoC;
        public Map<Integer, Double> tcPerReadCount;

        public TtoCData() {
            tcPerReadCount = new HashMap<>();
        }

        public double getTtoCRate() {
            return totalT > 0 ? totalTtoC/totalT : 0;
        }

        public double getNoTcReadCount() {
            if (!tcPerReadCount.containsKey(0)) {
                return 0.d;
            }
            return tcPerReadCount.get(0);
        }

        public double getTcReadCount() {
            double total = 0.d;
            for (Integer i : tcPerReadCount.keySet()) {
                if (i == 0) {
                    continue;
                }
                total += tcPerReadCount.get(i);
            }
            return total;
        }
    }

    public static TtoCData calculateTtoCcount(ReferenceSequence ref, GenomicRegion reg, GenomicRegionStorage<AlignedReadsData> cit, Genomic genomic, int condition, boolean onlyOverlap) {
        TtoCData ttoCData = new TtoCData();
        cit.ei(ref, reg).filter(f -> {
            return reg.contains(GenomicRegionPosition.FivePrime.position(f));
        }).filter(f-> {
            return !onlyOverlap || f.getData().getGeometryOverlap(0) >= 1;
        }).filter(f -> {
            for (int d = 0; d < f.getData().getDistinctSequences(); d++) {
                if (f.getData().getCount(d, condition, ReadCountMode.Weight) > 0) {
                    return true;
                }
            }
            return false;
        }).forEachRemaining(r -> {
            AlignedReadsData dat = r.getData();
            double totalCondCount = dat.getTotalCountForCondition(condition, ReadCountMode.Weight);
            double tcount = StringUtils.countChar(genomic.getSequence(r).subSequence(dat.getGeometryBeforeOverlap(0), dat.getGeometryBeforeOverlap(0) + dat.getGeometryOverlap(0)), 'T');
            if (!onlyOverlap) {
                tcount = StringUtils.countChar(genomic.getSequence(r), 'T');
            }
            ttoCData.totalT += (tcount * totalCondCount);
            if (tcount == 0) {
                return;
            }
            for (int d = 0; d < dat.getDistinctSequences(); d++) {
                AlignedReadsVariation[] variations = dat.getVariations(d);
                Set<Integer> pos = new HashSet<>();
                double readTcCount = dat.getCount(d, condition, ReadCountMode.Weight);
                int readTcNum = 0;
                for (int v = 0; v < variations.length; v++) {
                    if (!variations[v].isMismatch()) {
                        continue;
                    }
//                    if (snps != null && snps.contains(r.getReference().getName() + ":" + r.getRegion().map(variations[v].getPosition()))) {
//                        continue;
//                    }
                    if (variations[v].getReferenceSequence().length() > 1) {
                        System.err.println(variations[v].getReferenceSequence());
                    }
                    if (!variations[v].isFromSecondRead() && variations[v].getReferenceSequence().charAt(0) == 'T' && variations[v].getReadSequence().charAt(0) == 'C') {
                        if (pos.contains(variations[v].getPosition())) {
                            ttoCData.totalTtoC += dat.getCount(d, condition, ReadCountMode.Weight);
                            readTcNum++;
                            pos.remove(variations[v].getPosition());
                        } else if (!onlyOverlap && !dat.isPositionInOverlap(d, variations[v].getPosition())) {
                            ttoCData.totalTtoC += dat.getCount(d, condition, ReadCountMode.Weight);
                            readTcNum++;
                        } else {
                            pos.add(variations[v].getPosition());
                        }
                    }
                    if (variations[v].isFromSecondRead() && variations[v].getReferenceSequence().charAt(0) == 'A' && variations[v].getReadSequence().charAt(0) == 'G') {
                        if (pos.contains(variations[v].getPosition())) {
                            ttoCData.totalTtoC += dat.getCount(d, condition, ReadCountMode.Weight);
                            readTcNum++;
                            pos.remove(variations[v].getPosition());
                        } else if (!onlyOverlap && !dat.isPositionInOverlap(d, variations[v].getPosition())) {
                            ttoCData.totalTtoC += dat.getCount(d, condition, ReadCountMode.Weight);
                            readTcNum++;
                        } else {
                            pos.add(variations[v].getPosition());
                        }
                    }
                }
                double current = ttoCData.tcPerReadCount.getOrDefault(readTcNum, 0.0d);
                ttoCData.tcPerReadCount.put(readTcNum, current + readTcCount);
            }
        });
        return ttoCData;
    }

    // This function is used for the iTiSS-paper only and transforms data formats of ADAPT-CAGE, TSRFinder, TSSPredator into a uniform one
    // p-Values are switched (i.e. 1-pValue) to unify the meaning of scores (high totalDelta = TiSS)
    // cRNA-seq data is normalized for all three metrices (DENSITY/KINETIC/PEAK)
    //    This means, that the values are converted into 0-1 range, where 0 is the lowest value and 1 the highest value for the respective dataset
    //
    public static void transformData(String[] inFiles, String outFile, int window, boolean isPVal, int skip, int refCol, int strandCol, int tssCol, int valCol, boolean normalize, boolean twoVal) throws IOException {
        Map<ReferenceSequence, List<MutablePair<Integer, Double>>> map = new HashMap<>();
        for (String inFile : inFiles) {
            Map<ReferenceSequence, List<MutablePair<Integer, Double>>> mapTmp = new HashMap<>();
            EI.lines(inFile).skip(skip).forEachRemaining(l -> {
                String[] split = StringUtils.split(l, "\t");
                ReferenceSequence ref = Chromosome.obtain(split[refCol]);
                if (strandCol >= 0) {
                    ref = Chromosome.obtain(split[refCol] + split[strandCol]);
                }
                int pos = Integer.parseInt(split[tssCol]);
                double val;
                if (twoVal && split.length > 3) {
                    double v1 = Double.parseDouble(split[valCol]);
                    double v2 = Double.parseDouble(split[valCol+1]);
                    val = v1 > v2 ? v1 : v2;
                } else if (isPVal) {
                    val = 1 - Double.parseDouble(split[valCol]);
                } else {
                    if (split[valCol].contains(">")) {
                        val = 100d; // This needs to be done for TSSPredator data...
                    } else {
                        val = Double.parseDouble(split[valCol]);
                    }
                }
                mapTmp.computeIfAbsent(ref, absent -> new ArrayList<>()).add(new MutablePair<>(pos, val));
            });
            if (normalize) normalize(mapTmp);
            for (ReferenceSequence ref : mapTmp.keySet()) {
                map.computeIfAbsent(ref, absent -> new ArrayList<>()).addAll(mapTmp.get(ref));
            }
        }
        LineWriter writer = new LineOrientedFile(outFile).write();
        writer.writeLine("Ref\tTiSS\tValue");
        EI.wrap(map.keySet()).forEachRemaining(ref -> {
            List<MutablePair<Integer, Double>> lst = map.get(ref);
            lst.sort(Comparator.comparingInt(a -> a.Item1));

            int lastPos = lst.get(0).Item1;
            MutablePair<Integer, Double> lastHighest = lst.get(0);
            for (int i = 1; i < lst.size(); i++) {
                MutablePair<Integer, Double> current = lst.get(i);
                if (current.Item1-lastPos >= window) {
                    try {
                        writer.writeLine(ref + "\t" + lastHighest.Item1 + "\t" + lastHighest.Item2);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    lastHighest = current;
                } else {
                    lastHighest = current.Item2 > lastHighest.Item2 ? current : lastHighest;
                }
                lastPos = current.Item1;
            }
            try {
                writer.writeLine(ref + "\t" + lastHighest.Item1 + "\t" + lastHighest.Item2);
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        writer.close();
    }
}
