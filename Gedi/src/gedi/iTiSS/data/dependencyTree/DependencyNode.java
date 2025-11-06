package gedi.iTiSS.data.dependencyTree;

import java.util.ArrayList;
import java.util.List;

public class DependencyNode<T> {
    private List<DependencyNode<T>> inNodes;
    private List<DependencyNode<T>> outNodes;
    private int tsrThresh;
    private int tissThresh;
    private int id;
    private T data;
    private boolean checked;

    public DependencyNode(int id, List<DependencyNode<T>> inNodes, int tsrThresh, int tissThresh) {
        this.id = id;
        this.inNodes = inNodes;
        this.tsrThresh = tsrThresh;
        this.tissThresh = tissThresh;
    }

    public void addOutput(DependencyNode<T> out) {
        if (this.outNodes == null) {
            this.outNodes = new ArrayList<>();
        }
        this.outNodes.add(out);
    }

    public DependencyNode<T> getNodeWithId(int id) {
        if (this.id == id) {
            return this;
        }
        if (outNodes == null) {
            return null;
        }
        for (DependencyNode<T> node : outNodes) {
            DependencyNode<T> ret = node.getNodeWithId(id);
            if (ret != null) {
                return ret;
            }
        }
        return null;
    }

    public int getId() {
        return this.id;
    }

    public List<DependencyNode<T>> getOutNodes() {
        return outNodes;
    }

    public T getData() {
        return data;
    }

    public void setData(T data) {
        this.data = data;
    }

    public boolean isChecked() {
        return checked;
    }

    public void setChecked(boolean checked) {
        this.checked = checked;
    }

    public List<DependencyNode<T>> getInNodes() {
        return inNodes;
    }

    public int getTsrThresh() {
        return tsrThresh;
    }

    public int getTissThresh() {
        return tissThresh;
    }
}
