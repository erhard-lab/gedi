package gedi.iTiSS.merger2;

import gedi.core.reference.ReferenceSequence;

public class TissFileEntry {
    private ReferenceSequence reference;
    private int position;
    private Double readcount;
    private Double value;
    private int originId;

    public TissFileEntry(ReferenceSequence reference, int position, double readcount, double value, int originId) {
        this.reference = reference;
        this.position = position;
        this.readcount = readcount;
        this.value = value;
        this.originId = originId;
    }

    public void setOriginId(int originId) {
        this.originId = originId;
    }

    public ReferenceSequence getReference() {
        return reference;
    }

    public int getPosition() {
        return position;
    }

    public Double getReadcount() {
        return readcount;
    }

    public Double getValue() {
        return value;
    }

    public int getOriginId() {
        return originId;
    }
}
