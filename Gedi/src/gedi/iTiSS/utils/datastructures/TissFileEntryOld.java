package gedi.iTiSS.utils.datastructures;

import gedi.core.reference.ReferenceSequence;

public class TissFileEntryOld<T> {
    private ReferenceSequence reference;
    private int position;
    private T data;

    public TissFileEntryOld(ReferenceSequence reference, int position, T data) {
        this.reference = reference;
        this.position = position;
        this.data = data;
    }

    public ReferenceSequence getReference() {
        return reference;
    }

    public int getPosition() {
        return position;
    }

    public T getData() {
        return data;
    }

    public TissFileEntryOld<T> toOppositeStrand() {
        return new TissFileEntryOld<>(getReference().toOppositeStrand(), getPosition(), getData());
    }

    @Override
    public boolean equals(Object obj) {
        if (!(obj instanceof TissFileEntryOld)) {
            return false;
        }
        TissFileEntryOld other = (TissFileEntryOld) obj;
        return other.getPosition() == getPosition() &&
                other.getReference().equals(getReference());
    }

    @Override
    public int hashCode() {
        return (reference.getName().hashCode() << 28) | (reference.getStrand().ordinal() << 26) | position;
    }
}
