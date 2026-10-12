package gedi.util.datastructure.array;


import java.util.Arrays;
import java.util.EmptyStackException;

public class IntStack {
    private int[] data;
    private int size;

    public IntStack(int initialCapacity) {
        if (initialCapacity < 0) {
            throw new IllegalArgumentException("Negative initial capacity: " + initialCapacity);
        }
        this.data = new int[initialCapacity];
        this.size = 0;
    }

    public IntStack() {
        this(16);
    }

    public void push(int val) {
        if (size == data.length) {
            data = Arrays.copyOf(data, Math.max(data.length * 2, 4));
        }
        data[size++] = val;
    }

    public int pop() {
        if (size == 0) {
            throw new EmptyStackException();
        }
        return data[--size];
    }
    
    public int peek() {
        if (size == 0) {
            throw new EmptyStackException();
        }
        return data[size - 1];
    }

    public int size() {
        return size;
    }
}