package gedi.iTiSS.utils.loader;

import gedi.core.data.annotation.ScoreProvider;
import gedi.core.reference.ReferenceSequence;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class TsrData implements BinarySerializable, ScoreProvider {
    private int score;
    private int maxTissPos;
    private Set<Integer> calledBy;
    private Set<Integer> allTissPos;

    public TsrData(int score, int maxTissPos, Set<Integer> calledBy) {
        this.score = score;
        this.maxTissPos = maxTissPos;
        this.calledBy = calledBy;
    }

    public void setScore(int score) {
        this.score = score;
    }

    public void setMaxTissPos(int maxTissPos) {
        this.maxTissPos = maxTissPos;
    }

    public void setCalledBy(Set<Integer> calledBy) {
        this.calledBy = calledBy;
    }

    public void setAllTissPos(Set<Integer> allTissPos) {
        this.allTissPos = allTissPos;
    }

    public TsrData() {
        calledBy = new HashSet<>();
        allTissPos = new HashSet<>();
    }

    public Set<Integer> getAllTissPos() {
        return allTissPos;
    }

    public int getMaxTissPos() {
        return maxTissPos;
    }

    public Set<Integer> getCalledBy() {
        return calledBy;
    }

    public int getMostUpstreamTissPosition(ReferenceSequence ref) {
        IntStream allTissStream = allTissPos.stream().mapToInt(Integer::intValue);
        OptionalInt tissPos =  ref.isPlus() ? allTissStream.min() : allTissStream.max();
        if (!tissPos.isPresent()) throw new IllegalStateException("No TiSS in the set");
        return tissPos.getAsInt();
    }

    @Override
    public double getScore() {
        return score;
    }

    @Override
    public void serialize(BinaryWriter out) throws IOException {
        out.putInt(score);
        out.putInt(maxTissPos);
    }

    @Override
    public void deserialize(BinaryReader in) throws IOException {
        score = in.getInt();
        maxTissPos = in.getInt();
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        return (prime + score) * prime + maxTissPos;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        TsrData td = (TsrData) obj;
        return score == td.score && maxTissPos == td.maxTissPos;
    }

    @Override
    public String toString() {
        return String.format("Max TiSS position: %d, Score: %d", maxTissPos, score);
    }
}
