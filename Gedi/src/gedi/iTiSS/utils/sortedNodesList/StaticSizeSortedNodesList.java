package gedi.iTiSS.utils.sortedNodesList;

import java.util.*;

public class StaticSizeSortedNodesList<T> {
    private Node<T> first;
    private Comparator<T> comparator;

    /**
     * Creates a sorted linked list static in size. It is sorted based on the {@code comparator}.
     * The special feature of this list is its ability to getValueAtIndex references to single {@code Node}s as well as
     * to retrieve its values like usual.
     * {@code NULL} values are not allowed.
     * @param collection Will be transformed into a {@code StaticSizeSortedNodesList}
     * @param comparator Used for sorting
     */
    public StaticSizeSortedNodesList(Collection<T> collection, Comparator<T> comparator) {
        this.comparator = comparator;
        collection.forEach(this::insert);
    }

    private void insert(T value) {
        if (value == null) {
            throw new IllegalArgumentException("NULL values are not supported for StaticSortedNodesList");
        }
        Node<T> node = new Node<>(value);
        if (first == null) {
            first = node;
        } else if (comparator.compare(first.getValue(), value) >= 0) {
            node.setNext(first);
            first = node;
        } else {
            Node<T> current = first;
            while (current.hasNext()) {
                if (comparator.compare(current.getNext().getValue(), value) >= 0) {
                    node.setNext(current.getNext());
                    current.setNext(node);
                    return;
                }
            }
            current.setNext(node);
        }
    }

    private boolean delete(T value) {
        if (value == null) {
            throw new IllegalArgumentException("NULL values are not supported for StaticSortedNodesList");
        }
        if (first == null) {
            throw new IllegalArgumentException("Empty list. Should never happen.");
        }
        if (comparator.compare(first.getValue(), value) == 0) {
            first = first.getNext();
            return true;
        }
        Node<T> current = first;
        while (current.hasNext()) {
            if (comparator.compare(current.getNext().getValue(), value) == 0) {
                current.setNext(current.getNext().getNext());
                return true;
            }
        }
        return false;
    }

    public boolean insertAndDelete(T valToInsert, T valToDelete) {
        if (delete(valToDelete)) {
            insert(valToInsert);
            return true;
        }
        return false;
    }

    public ImmutableNode<T> getNodeReference(int index) {
        Node<T> returnNode = first;
        while (index != 0) {
            returnNode = returnNode.getNext();
            index--;
        }
        return returnNode;
    }
}
