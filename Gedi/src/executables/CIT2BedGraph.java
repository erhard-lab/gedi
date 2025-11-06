package executables;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.iTiSS.utils.GenomicUtils;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.iTiSS.utils.WriterUtils;
import gedi.iTiSS.utils.datastructures.SparseNumericArray;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.userInteraction.progress.ConsoleProgress;

import java.io.IOException;
import java.util.Comparator;

public class CIT2BedGraph {
    public static void main(String[] args) throws IOException {
        if (args.length < 4 || !args[args.length-2].endsWith(".cit")) {
            usage();
            return;
        }
        boolean showProgress = false;
        ConversionType covType = ConversionType.Coverage;
        Strandness strandness = Strandness.Sense;
        Genomic genomic = null;
        boolean onlyStandardRefs = false;
        boolean totalize = false;
        int i = 0;
        for (; i < args.length; i++) {
            if (args[i].equals("-p")) {
                showProgress = true;
            } else if (args[i].equals("-t")) {
                i++;
                ConversionType.valueOf(args[i]);
            } else if (args[i].equals("-s")) {
                i++;
                strandness = Strandness.valueOf(args[i]);
            } else if (args[i].equals("-g")) {
                i++;
                genomic = Genomic.get(args[i]);
            } else if (args[i].equals("-stdChr")) {
                onlyStandardRefs = true;
            } else if (args[i].equals("-total")) {
                totalize = true;
            } else {
                break;
            }
        }
        if (i != args.length-2) {
            usage();
            return;
        }
        GenomicRegionStorage<AlignedReadsData> reads = new CenteredDiskIntervalTreeStorage<>(args[i]);
        i++;
        String[] sampleNames = reads.getMetaDataConditions();

        writeCountsToBedGraph(args[i], sampleNames, reads, genomic, strandness, true, showProgress, onlyStandardRefs, totalize);
        writeCountsToBedGraph(args[i], sampleNames, reads, genomic, strandness, false, showProgress, onlyStandardRefs, totalize);

        final Genomic fGenomic = genomic;
        LineWriter chrSizeWriter = new LineOrientedFile(args[i] + "chr.sizes").write();
        reads.iterateReferenceSequences()
                .iff(onlyStandardRefs, it -> it.filter(GenomicUtils::isStandardChromosome))
                .map(ReferenceSequence::toStrandIndependent).unique(true)
                .map(ref -> GenomicUtils.toChromosomeAddedName(ref) + "\t" + fGenomic.getLength(ref.getName()))
                .print(chrSizeWriter);
        chrSizeWriter.close();
    }

    public static void writeCountsToBedGraph(String prefix, String[] sampleNames, GenomicRegionStorage<AlignedReadsData> reads,
                                             Genomic fGenomic, Strandness fStrandness, boolean isPlus, boolean showProgress,
                                             boolean onlyStandardRefs, boolean totalize) {
        LineWriter[] plus_writers = totalize ? new LineWriter[] {new LineOrientedFile(isPlus ? prefix + "_plus.bedGraph" : prefix + "_minus.bedGraph").write()} : WriterUtils.createWriters(prefix, sampleNames, isPlus ? "_plus.bedGraph" : "_minus.bedGraph");
        int count = reads.iterateReferenceSequences().filter(ref -> ref.isPlus() == isPlus).countInt();
        reads.iterateReferenceSequences()
                .iff(onlyStandardRefs, it -> it.filter(GenomicUtils::isStandardChromosome))
                .sort(Comparator.comparing(GenomicUtils::toChromosomeAddedName))
                .filter(ref -> ref.isPlus() == isPlus)
                .iff(showProgress, it -> it.progress(new ConsoleProgress(System.err), count, ReferenceSequence::toPlusMinusString))
                .forEachRemaining(ref -> {
            SparseNumericArray<Float>[] ary = new SparseNumericArray[sampleNames.length];
            for (int cond = 0; cond < sampleNames.length; cond++) {
                ary[cond] = new SparseNumericArray<>(fGenomic.getLength(ref.toPlusMinusString()), 0.0f);
            }
            TiSSUtils.extractFivePrimeCountsFromSingleFileSparseFull(ary, reads, ref, fStrandness);
            ary = totalize ? totalize(ary) : ary;
            for (int cond = 0; cond < ary.length; cond++) {
                for (int pos : EI.wrap(ary[cond].getNonZeroIndices()).sort(Comparator.comparingInt(comp -> comp)).loop()) {
                    plus_writers[cond].writeLine2(GenomicUtils.toChromosomeAddedName(ref) + "\t" + (pos+1) + "\t" + (pos+2) + "\t" + ary[cond].get(pos)); // 0-index to 1-index
                }
            }
        });
        WriterUtils.closeWriters(plus_writers);
    }

    private static SparseNumericArray<Float>[] totalize(SparseNumericArray<Float>[] ary) {
        SparseNumericArray<Float> total = new SparseNumericArray<>(ary[0].length(), 0.0f);
        for (int cond = 0; cond < ary.length; cond++) {
            for (int pos : EI.wrap(ary[cond].getNonZeroIndices()).loop()) {
                total.set(pos, ary[cond].get(pos) + total.get(pos));
            }
        }
        return new SparseNumericArray[] {total};
    }

    public static void usage() {
        System.out.println("CIT2BedGraph [Options] <in.cit> <outPrefix>");
        System.out.println("\tOptions:");
        System.out.println("\t\t-p\t\t\tShow progress");
        System.out.println("\t\t-s\t\t\tStrandness [Sense, Antisense]");
        System.out.println("\t\t-g\t\t\tGenomic");
        System.out.println("\t\t-stdChr\t\tUse only standard chromosomes (1,2,...,X,Y)");
        System.out.println("\t\t-t\t\t\tConversion type [FivePrime, ThreePrime, Coverage]");
        System.out.println("\t\t-total\t\tTotalize all conditions into one");
    }

    public enum ConversionType { FivePrime, ThreePrime, Coverage }
}
