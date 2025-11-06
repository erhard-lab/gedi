package gedi.iTiSS.utils.sortedNodesList;

public class Node<T> extends ImmutableNode<T> {
    private Node<T> next;

    public Node(T value) {
        this.value = value;
    }

    public void setNext(Node<T> next) {
        this.next = next;
    }

    public Node<T> getNext() {
        return next;
    }

    public boolean hasNext() {
        return next != null;
    }
}
