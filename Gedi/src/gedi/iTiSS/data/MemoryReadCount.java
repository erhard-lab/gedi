package gedi.iTiSS.data;

import gedi.util.datastructure.array.NumericArray;

public class MemoryReadCount {
    private NumericArray[] readCounts;
    private int finishedAccesses;
    private int maxAccessCount;
    private boolean safeToDelete;
    private boolean forceDeleted;
    public int currentAccessCount;
    private boolean multi;
    private boolean inMemory;
    private boolean loading;
    private long byteSize;

    public MemoryReadCount(int maxAccessCount, boolean multi) {
        this.maxAccessCount = maxAccessCount;
        this.multi = multi;
        finishedAccesses = 0;
    }

    public NumericArray[] startAccess(boolean multi) {
        if (!isInMemory()) {
            throw new IllegalArgumentException("Read counts not in memory. Wrongly accessed.");
        }
        currentAccessCount++;
        if (this.multi && !multi) {
            NumericArray total = readCounts[0].copy();
            for (int i = 1; i < readCounts.length; i++) {
                total.add(readCounts[i]);
            }
            return new NumericArray[] {total};
        }
        return readCounts;
    }

    public void finishAccess() {
        finishedAccesses++;
        currentAccessCount--;
        if (finishedAccesses >= maxAccessCount) {
            safeToDelete = true;
            removeFromMemory();
        }
    }

    public boolean forceDelete() {
        if (!isInMemory()) {
            return false;
        }
        if (currentAccessCount == 0) {
            removeFromMemory();
            forceDeleted = true;
            return true;
        }
        return false;
    }

    private void removeFromMemory() {
        inMemory = false;
        for (int i = 0; i < readCounts.length; i++) {
            readCounts[i] = null;
        }
        readCounts = null;
        byteSize = 0;
    }

    public void setReadCount(NumericArray[] readCounts) {
        this.readCounts = readCounts;
        this.byteSize = DataWrapper.SINGLE_VALUE_BYTE_SIZE * readCounts[0].length() * readCounts.length + readCounts.length;
        inMemory = true;
        loading = false;
    }

    public void setLoading(boolean loading) {
        this.loading = loading;
    }

//    public boolean isInMemory() {
//        if (readCounts == null) {
//            return false;
//        }
//        for (NumericArray ary : readCounts) {
//            if (ary == null) {
//                return false;
//            }
//        }
//        return true;
//    }

    public boolean isInMemory() {
        return inMemory;
    }

    public boolean isSafeToDelete() {
        return safeToDelete;
    }

    public boolean isMulti() {
        return multi;
    }

    public boolean isLoading() {
        return loading;
    }

    public long getByteSize() {
        return byteSize;
    }
}
