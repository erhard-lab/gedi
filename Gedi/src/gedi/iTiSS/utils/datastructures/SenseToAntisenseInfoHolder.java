package gedi.iTiSS.utils.datastructures;

public class SenseToAntisenseInfoHolder {
    private TissFileEntryOld<Double> senseTiss;
    private TissFileEntryOld<Double> antisenseTiss;
    private String sequence;

    public SenseToAntisenseInfoHolder(TissFileEntryOld<Double> senseTiss, TissFileEntryOld<Double> antisenseTiss, String sequence) {
        this.senseTiss = senseTiss;
        this.antisenseTiss = antisenseTiss;
        this.sequence = sequence;
    }

    // Negative = upstream, Positive = downstream
    public int getDistance() {
        if (senseTiss.getReference().isMinus()) {
            return senseTiss.getPosition() - antisenseTiss.getPosition();
        } else if (senseTiss.getReference().isPlus()) {
            return antisenseTiss.getPosition() - senseTiss.getPosition();
        } else {
            return Math.abs(antisenseTiss.getPosition() - senseTiss.getPosition());
        }
    }

    public String getSequence() {
        return sequence;
    }

    public String getFastaId() {
        return senseTiss.getReference().toPlusMinusString() + "_" +
                senseTiss.getPosition() + "_" + senseTiss.getData() + " | " +
                antisenseTiss.getReference().toPlusMinusString() + "_" +
                antisenseTiss.getPosition() + "_" + antisenseTiss.getData() + " | " +
                getDistance();
    }

    public double getSenseEnrichment() {
        return senseTiss.getData();
    }

    public double getAntisenseEnrichment() {
        return antisenseTiss.getData();
    }

    public static String getStatsHeader() {
        return "SenseRef\tAntisenseRef\tSensePos\tAntisensePos\tSenseEnrichment\tAntisenseEnrichment\tDistance";
    }

    @Override
    public String toString() {
        return senseTiss.getReference().toPlusMinusString() + "\t" +
                antisenseTiss.getReference().toPlusMinusString() + "\t" +
                senseTiss.getPosition() + "\t" +
                antisenseTiss.getPosition() + "\t" +
                senseTiss.getData() + "\t" +
                antisenseTiss.getData() + "\t" +
                getDistance();
    }
}
