package gedi.iTiSS.data.dependencyTree;

public class DependencyNodeMeta {
    private int[] input;
    private int output;
    private int tsrThresh;
    private boolean added;
    private int tissThresh;

    public DependencyNodeMeta(int[] input, int output, int tsrThresh, int tissThresh) {
        this.input = input;
        this.output = output;
        this.tsrThresh = tsrThresh;
        this.tissThresh = tissThresh;
    }

    public int[] getInput() {
        return input;
    }

    public int getOutput() {
        return output;
    }

    public int getTsrThresh() {
        return tsrThresh;
    }

    public int getTissThresh() {
        return tissThresh;
    }

    public void add() {
        added = true;
    }

    public boolean isAdded() {
        return added;
    }
}
