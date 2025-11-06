package gedi.iTiSS.utils.sortedNodesList;

import java.util.*;

/**
 * Ascending ordered list
 * @param <T>
 */
public class StaticSizeSortedArrayList<T> {
    private List<T> list;
    private Comparator<T> comparator;
    private int size;


    public StaticSizeSortedArrayList(Collection<T> collection, Comparator<T> comparator) {
        this.comparator = comparator;
        list = new ArrayList<>(collection.size()*2);
        list.addAll(collection);
        list.sort(comparator);
        this.size = collection.size();
    }

    private void insert(T value) {
        if (value == null) {
            throw new IllegalArgumentException("NULL values are not supported for StaticSortedNodesList");
        }
        int index = Collections.binarySearch(list, value, comparator);
        if (index<0) {
            index = Math.abs(index)-1;
        }
        if (index == list.size()) {
            list.add(value);
        } else {
            list.add(index, value);
        }
    }

    private boolean delete(T value) {
        int index = Collections.binarySearch(list, value, comparator);
        if (index < 0) {
            return false;
        }
        list.remove(index);
        return true;
    }

    public boolean insertSortedAndDelete(T val2Insert, T val2Delete) {
        if (delete(val2Delete)) {
            insert(val2Insert);
            return true;
        }
        return false;
    }

    public T getValueAtIndex(int index) {
        return list.get(index);
    }

    /**
     * Returns the index of the next greater value than {@code value}.
     * Returns {@code list}.size() if the value is equal or greater than the largest value in {@code list}
     * @param value The value to check
     * @return Highest position of {@code value}
     */
    public int getNextHigherIndex(T value) {
        int index = Collections.binarySearch(list, value, comparator);
        if (index < 0) {
            index = Math.abs(index) - 1;
        }
        for (int i = index; i < list.size(); i++) {
            if (comparator.compare(list.get(i), value) > 0) {
                return i;
            }
        }
        return list.size();
    }

    /**
     * Returns the index of the next lower value than {@code value}.
     * Returns -1 if the value is equal or lower than the lowest value in {@code list}
     * @param value The value to check
     * @return Lowest position of {@code value}
     */
    public int getNextLowerIndex(T value) {
        int index = Collections.binarySearch(list, value, comparator);
        if (index < 0) {
            index = Math.abs(index) - 1;
        }
        if (index == list.size()) {
            return list.size()-1;
        }
        for (int i = index; i >= 0; i--) {
            if (comparator.compare(list.get(i), value) < 0) {
                return i;
            }
        }
        return -1;
    }

    public int getSize() {
        return size;
    }

    public List<T> toList() {
        return new ArrayList<>(list);
    }
}
