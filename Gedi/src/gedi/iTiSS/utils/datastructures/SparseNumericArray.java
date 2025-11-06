package gedi.iTiSS.utils.datastructures;

import gedi.util.datastructure.array.NumericArray;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

public class SparseNumericArray<T extends Number> {
    private Map<Integer, T> data;
    private T defaultValue;
    private int length;
    private Predicate<T> isZero;

    public SparseNumericArray(int length, T defaultValue) {
        this.length = length;
        this.data = new HashMap<>((int) (length*0.01), 1.0f);
        this.defaultValue = defaultValue;
        this.isZero = this::isZeroEntry;
    }

    public SparseNumericArray(int length, T defaultValue, Predicate<T> isZero) {
        this.length = length;
        this.data = new HashMap<>((int) (length * 0.1), 0.9f);
        this.defaultValue = defaultValue;
        this.isZero = isZero;
    }

    public void set(int index, T value) {
        if (index >= length()) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        if (isZero.test(value)) {
            return;
        }
        data.put(index, value);
    }

    private boolean isZeroEntry(T value) {
        return value.floatValue() == 0.0f;
    }

    public Set<Integer> getNonZeroIndices() {
        return data.keySet();
    }

    public T get(int index) {
        if (index >= length()) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return data.getOrDefault(index, defaultValue);
    }

    public int length() {
        return length;
    }
}
