package gedi.iTiSS.data;

public class Data {
    private int[] lane;
    private boolean multi;

    public Data(int[] lane, boolean multi) {
        this.lane = lane;
        this.multi = multi;
    }

    public int[] getLane() {
        return lane;
    }

    public boolean isMulti() {
        return multi;
    }
}
