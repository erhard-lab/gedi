package gedi.iTiSS.data;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.iTiSS.utils.ReadType;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.riboseq.utils.RiboUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;

import java.util.*;

/**
 * Container holding and accessing the raw data.
 *
 * This class contains the data and serves as the main point of interacting with the data.
 * The raw data should never be accessed directly anywhere else than here.
 * This class provides functions to access the raw data indirectly.
 *
 * @author Christopher Juerges
 * @version 1.0
 */
public class DataWrapper {
    /**
     * The raw data. Basically the cit-files concatenated in a list ordered in the way the user entered them
     * in the arguments
     */
    public final static long SINGLE_VALUE_BYTE_SIZE = 6;
    public final static long GB = 1024*1024*1024;

    private List<GenomicRegionStorage<AlignedReadsData>> rawData;
    private Strandness strandness;
    private ReadType readType;
    private Set<ReferenceSequence> testChr;

    private Map<Set<Integer>, Map<ReferenceSequence, MemoryReadCount>> memoryMap;

    private List<MemoryReadCount> memoryReadCountAccessionOrder = new ArrayList<>();

    private long vmMemory = Runtime.getRuntime().maxMemory();
    private long currentlyUsedVmMemory = GB;

    private Set<ReferenceSequence> loadedChromosomes;

    public DataWrapper(List<GenomicRegionStorage<AlignedReadsData>> rawData, Strandness strandness, ReadType readType, Set<ReferenceSequence> testChr) {
        this.rawData = rawData;
        this.strandness = strandness;
        this.readType = readType;
        this.testChr = testChr;
    }

    public DataWrapper(List<GenomicRegionStorage<AlignedReadsData>> rawData, Strandness strandness) {
        this(rawData, strandness, ReadType.FIVE_PRIME, null);
    }

    public void initData(Genomic genomic, List<Data> lanes) {
        memoryMap = new HashMap<>();
        Map<Set<Integer>, Integer> accessCounts = calculateAccessCounts(lanes);
        loadedChromosomes = new HashSet<>();

        for (Set<Integer> key : accessCounts.keySet()) {
            Map<ReferenceSequence, MemoryReadCount> refMemoryMap = new HashMap<>();
            boolean multi = hasMulti(lanes, key);
            genomic.iterateReferenceSequences().forEachRemaining(r -> {
                if (testChr != null && !testChr.contains(r)) {
                    return;
                }
                loadedChromosomes.add(r);
                refMemoryMap.put(r, new MemoryReadCount(accessCounts.get(key), multi));
            });
            memoryMap.put(key, refMemoryMap);
        }
    }

    private boolean hasMulti(List<Data> lanes, Set<Integer> data) {
        for (Data d : lanes) {
            if (EI.wrap(d.getLane()).set().equals(data) && d.isMulti()) {
                return true;
            }
        }
        return false;
    }

    private Map<Set<Integer>, Integer> calculateAccessCounts(List<Data> lanes) {
        Map<Set<Integer>, Integer> accessCounts = new HashMap<>();
        for (Data l : lanes) {
            Set<Integer> lst = EI.wrap(l.getLane()).set();
            boolean known = accessCounts.containsKey(lst);
            if (known) {
                accessCounts.put(lst, accessCounts.get(lst) + 1);
            } else {
                accessCounts.put(lst, 1);
            }
        }
        return accessCounts;
    }

    private synchronized long getFreeMemoryIn() {
        return (vmMemory-Runtime.getRuntime().totalMemory())+Runtime.getRuntime().freeMemory();
    }

    private boolean freeMemory(long neededMemory) {
//        System.err.println("[" + Thread.currentThread().getName() + "] " + "Starting force-delete. list size: " + memoryReadCountAccessionOrder.size());
        for (int i = memoryReadCountAccessionOrder.size()-1; i >= 0; i--) {
            if (memoryReadCountAccessionOrder.get(i).getByteSize() > 0) {
//                System.err.println("[" + Thread.currentThread().getName() + "] " + i + ": access: " + memoryReadCountAccessionOrder.get(i).currentAccessCount + ", size: " + (memoryReadCountAccessionOrder.get(i).getByteSize() / (double) GB) + ", memory: " + (memoryReadCountAccessionOrder.get(i).isInMemory() ? "Y" : "N"));
            }
            long byteSize = memoryReadCountAccessionOrder.get(i).getByteSize();
            if (byteSize > 0 && memoryReadCountAccessionOrder.get(i).forceDelete()) {
                currentlyUsedVmMemory -= byteSize;
            }
//            if (currentlyUsedVmMemory + neededMemory < vmMemory && getFreeMemoryIn() > 0) {
            if (neededMemory < getFreeMemoryIn() && getFreeMemoryIn() > vmMemory/2) {
                return true;
            }
        }
        return false;
    }

    private void printMemoryUsage(long neededMemory) {
        System.err.println(String.format("[" + Thread.currentThread().getName() + "] " + "**** c: %.2f Gb, m: %.2f Gb, n: %.2f Gb, f: %.2f Gb, t: %.2f Gb",
                ((double)currentlyUsedVmMemory/(double)GB),
                ((double)vmMemory/(double)GB),
                ((double)neededMemory/(double)GB),
                ((double)getFreeMemoryIn()/(double)GB),
                ((double)Runtime.getRuntime().totalMemory()/(double)GB)));
    }

    private synchronized boolean checkAndFreeMemory(long neededMemory) {
        neededMemory *= 1.1;
//        printMemoryUsage(neededMemory);
        if ((double)neededMemory> vmMemory) {
            throw new OutOfMemoryError("Your data is too big for your VMs memory capacity. Either run JavaVM with more memory or totalize your data before.");
        }
        if (currentlyUsedVmMemory+neededMemory > vmMemory) {
//            System.err.println("[" + Thread.currentThread().getName() + "] " + "*** Freeing memory ***");
        } else {
            return true;
        }
        if (freeMemory(neededMemory)) {
//            System.err.println("[" + Thread.currentThread().getName() + "] " + "Memory successfully freed.");
            return true;
        }
        return false;
    }

    private synchronized void addUsedMemory(long byteSizeNeeded) {
        currentlyUsedVmMemory += byteSizeNeeded;
    }

    public synchronized NumericArray[] startAccessingData(Data data, ReferenceSequence ref, int refLength) {
        Set<Integer> laneList = EI.wrap(data.getLane()).set();
        MemoryReadCount memoryReadCount = memoryMap.get(laneList).get(ref);
        if (!memoryReadCount.isInMemory()) {
            long neededMemory = SINGLE_VALUE_BYTE_SIZE * refLength * data.getLane().length;
//            System.err.println(Thread.currentThread().getName() + " askes for memory");
            boolean freedSuccessfully = checkAndFreeMemory(neededMemory);
            if (!freedSuccessfully) {
                return null;
            }
            if (memoryReadCount.isMulti()) {
                memoryReadCount.setReadCount(loadMultiReadCountToMemory(data.getLane(), ref, refLength));
            } else {
                memoryReadCount.setReadCount(new NumericArray[] {totalizeReadCounts(loadMultiReadCountToMemory(data.getLane(), ref, refLength))});
            }
            memoryReadCountAccessionOrder.add(memoryReadCount);
            addUsedMemory(memoryReadCount.getByteSize());
        }
        return memoryReadCount.startAccess(data.isMulti());
    }

    public synchronized NumericArray[] startAccessingData2(Data data, ReferenceSequence ref, int refLength) {
        Set<Integer> laneList = EI.wrap(data.getLane()).set();
        MemoryReadCount memoryReadCount = memoryMap.get(laneList).get(ref);
        if (!memoryReadCount.isInMemory()) {
            if (readType == ReadType.DENSITY) {
                if (memoryReadCount.isMulti()) {
                    throw new RuntimeException();
                } else {
                    memoryReadCount.setReadCount(new NumericArray[] {TiSSUtils.extractReadDensities(rawData, data.getLane(), ref, refLength, strandness)});
                }
            } else {
                if (memoryReadCount.isMulti()) {
                    memoryReadCount.setReadCount(loadMultiReadCountToMemory(data.getLane(), ref, refLength));
                } else {
                    memoryReadCount.setReadCount(new NumericArray[]{TiSSUtils.extractCounts(rawData, data.getLane(), ref, refLength, strandness, readType)});
                }
            }
            memoryReadCountAccessionOrder.add(memoryReadCount);
            addUsedMemory(memoryReadCount.getByteSize());
        }
        return memoryReadCount.startAccess(data.isMulti());
    }

    public synchronized void finishAccessingData(Data lane, ReferenceSequence ref) {
        MemoryReadCount memoryReadCount = memoryMap.get(EI.wrap(lane.getLane()).set()).get(ref);
        long byteSize = memoryReadCount.getByteSize();
        memoryReadCount.finishAccess();
        if (memoryReadCount.getByteSize() == 0) {
            printMemoryUsage(-byteSize);
            currentlyUsedVmMemory -= byteSize;
        }
    }

    private NumericArray[] loadMultiReadCountToMemory(int[] lane, ReferenceSequence ref, int refLength) {
        CitAccessInfo citAccessInfo = getCitIndexAccessListNew(lane);
        NumericArray[] readCounts = new NumericArray[lane.length];
        int index = 0;
        for (int i = 0; i < citAccessInfo.getCitAccessNum(); i++) {
            NumericArray[] tmpReadCounts = loadMultiReadCountsToMemoryFromSingleFile(citAccessInfo.getLaneAccess(i), rawData.get(citAccessInfo.getCitAccess(i)), ref, refLength);
            for (NumericArray tmpReadCount : tmpReadCounts) {
                readCounts[index] = tmpReadCount;
                index++;
            }
        }
        return readCounts;
    }

    private NumericArray[] loadMultiReadCountsToMemoryFromSingleFile(List<Integer> citAccessIndeces, GenomicRegionStorage<AlignedReadsData> cit, ReferenceSequence ref, int refLength) {
        NumericArray[] readCounts = new NumericArray[citAccessIndeces.size()];
        for (int i = 0; i < readCounts.length; i++) {
            readCounts[i] = NumericArray.createMemory(refLength, NumericArray.NumericArrayType.Double);
        }
        final boolean switchFiveAndThreePrimeEnd = strandness == Strandness.Antisense && readType == ReadType.FIVE_PRIME || strandness == Strandness.Sense && readType == ReadType.THREE_PRIME;
        cit.ei(strandness.equals(Strandness.Antisense) ? ref.toOppositeStrand() : ref).forEachRemaining(r -> {
            int pos0 = switchFiveAndThreePrimeEnd ? GenomicRegionPosition.ThreePrime.position(r) : GenomicRegionPosition.FivePrime.position(r);
            int pos1 = switchFiveAndThreePrimeEnd ? GenomicRegionPosition.ThreePrime.position(r, 1) : GenomicRegionPosition.FivePrime.position(r, 1);
            NumericArray c0 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            NumericArray c1 = NumericArray.createMemory(r.getData().getNumConditions(), NumericArray.NumericArrayType.Double);
            for (int k = 0; k < r.getData().getDistinctSequences(); k++) {
                if (hasEndMismatch(r.getData(), k, r.getRegion().getTotalLength())) {
                    c1 = r.getData().addCountsForDistinct(k, c1, ReadCountMode.Weight);
                } else {
                    c0 = r.getData().addCountsForDistinct(k, c0, ReadCountMode.Weight);
                }
            }

            if (c0.length() > 0) {
                int index = 0;
                for (int i : citAccessIndeces) {
                    readCounts[index].setFloat(pos0, readCounts[index].getFloat(pos0) + c0.getFloat(i));
                    index++;
                }
            }
            if (c1.length() > 0) {
                int index = 0;
                for (int i : citAccessIndeces) {
                    readCounts[index].setFloat(pos1, readCounts[index].getFloat(pos1) + c1.getFloat(i));
                    index++;
                }
            }
        });
        return readCounts;
    }

    private boolean hasEndMismatch(AlignedReadsData ard, int distict, int readLength) {
        if (strandness.equals(Strandness.Antisense)) {
            return hasTailingMismatch(ard, distict, readLength);
        } else {
            return RiboUtils.hasLeadingMismatch(ard, distict);
        }
    }

    private boolean hasTailingMismatch(AlignedReadsData ard, int distinct, int readLength) {
        for (int i=0; i<ard.getVariationCount(distinct); i++) {
            if (ard.isMismatch(distinct, i) && ard.getMismatchPos(distinct, i)==readLength-1)
                return true;
        }
        return false;
    }

    private NumericArray totalizeReadCounts(NumericArray[] readCounts) {
        NumericArray totalizedReadCounts = readCounts[0].copy();
        for (int i = 1; i < readCounts.length; i++) {
            totalizedReadCounts.add(readCounts[i]);
        }
        return totalizedReadCounts;
    }

    private CitAccessInfo getCitIndexAccessListNew(int[] lane) {
        CitAccessInfo citAccessInfo = new CitAccessInfo();
        for (int l : lane) {
            int passedCondNum = 0;
            for (int i = 0; i < rawData.size(); i++) {
                int numCond = rawData.get(i).getMetaDataConditions().length;
                if (l < passedCondNum+numCond) {
                    citAccessInfo.add(i, l-passedCondNum);
                    break;
                }
                passedCondNum += numCond;
            }
        }
        return citAccessInfo;
    }

    public Set<ReferenceSequence> getLoadedChromosomes() {
        return new HashSet<>(loadedChromosomes);
    }
}
