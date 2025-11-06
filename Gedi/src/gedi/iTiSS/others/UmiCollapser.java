package gedi.iTiSS.others;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.*;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.userInteraction.progress.ConsoleProgress;

import java.io.IOException;
import java.util.*;

public class UmiCollapser {
    public static void collapsUmis(String path, String genomicName) throws IOException {
        Genomic genomic = Genomic.get(genomicName);
        CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> reads = new CenteredDiskIntervalTreeStorage<>(path);
        CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(path+".collapsed.cit", DefaultAlignedReadsData.class);

        ExtendedIterator<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> it = EI.empty();
        for (ReferenceSequence refSeq : genomic.iterateReferenceSequences().loop()) {
            it = it.chain(reads.ei(refSeq).progress(new ConsoleProgress(System.err), reads.ei(refSeq).countInt(), pr->"Processing " + pr.getReference().toPlusMinusString()).map(r -> {
                if (r.getRegion().getTotalLength() < 18) {
                    return null;
                }
                DefaultAlignedReadsData data = r.getData();
                return new ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>(r.getReference(), r.getRegion(), collapsUmiForRead(data));
            }).removeNulls());
        }

        out.fill(it);
    }

    private static DefaultAlignedReadsData collapsUmiForRead(DefaultAlignedReadsData data) {
        HashMap<String, List<Integer>> umi2Distincts = new HashMap<>();
        for (int i = 0; i < data.getDistinctSequences(); i++) {
            String umi = getSoftclip(data, i,false).getReadSequence().toString() + getSoftclip(data, i, true).getReadSequence().toString();
            if (umi2Distincts.keySet().contains(umi)) {
                umi2Distincts.get(umi).add(i);
                continue;
            }
            List<Integer> distincts = new ArrayList<>();
            umi2Distincts.put(umi, distincts);
            distincts.add(i);
        }

        AlignedReadsDataFactory factory = new AlignedReadsDataFactory(data.getNumConditions());
        factory.start();
        EI.wrap(umi2Distincts.keySet()).forEachRemaining(umi -> {
            factory.newDistinctSequence();
            int[] count = new int[data.getNumConditions()];
            int multiplicity = 0;
            for (int d : umi2Distincts.get(umi)) {
                for (int j = 0; j < count.length; j++) {
                    count[j] = data.getCount(d, j) + count[j] > 0 ? 1 : 0;
                }
                multiplicity += data.getMultiplicity(d);
            }
            factory.setCount(count);
            factory.setMultiplicity(multiplicity);
            int d = umi2Distincts.get(umi).get(0);
            for (int v = 0; v < data.getVariationCount(d); v++) {
                if (data.getVariation(d, v).isSoftclip()) {
                    factory.addVariation(data.getVariation(d, v));
                }
            }
            factory.setGeometry(data.getGeometryBeforeOverlap(d), data.getGeometryOverlap(d), data.getGeometryAfterOverlap(d));
        });
        return factory.create();
    }

    // TODO: Overlap changes not implemented yet!!!!!!!!!
    public static void cutUmi(String path, int umiLength, String genomicName) throws IOException {
        Genomic genomic = Genomic.get(genomicName);

        CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> reads = new CenteredDiskIntervalTreeStorage<>(path);
        CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(path+".umiCorrected.cit", DefaultAlignedReadsData.class);

        ExtendedIterator<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> it = EI.empty();
        for (ReferenceSequence refSeq : genomic.iterateReferenceSequences().loop()) {
            System.err.println(refSeq);
            it = it.chain(reads.ei(refSeq).progress(new ConsoleProgress(System.err), reads.ei(refSeq).countInt(), pr->"Processing " + pr.getReference().toPlusMinusString()).unfold(r -> {
                List<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> outList = new ArrayList<>();
                DefaultAlignedReadsData data = r.getData();
                AlignedReadsDataFactory factory = new AlignedReadsDataFactory(data.getNumConditions());
                factory.start();
                for (int i = 0; i < data.getDistinctSequences(); i++) {
                    GenomicRegion regNew = r.getRegion();
                    boolean regionChanged = false;
                    AlignedReadsSoftclip fPrime = getSoftclip(data, i, false);
                    AlignedReadsSoftclip tPrime = getSoftclip(data, i, true);
                    int fDelta = fPrime == null ? -umiLength : fPrime.getReadSequence().length() - umiLength;
                    int tDelta = tPrime == null ? -umiLength : tPrime.getReadSequence().length() - umiLength;
                    if (fDelta < 0) {
                        regNew = r.getReference().isPlus() ? regNew.extendFront(fDelta) : regNew.extendBack(fDelta);
                        regionChanged = true;
                    }
                    if (tDelta < 0) {
                        regNew = r.getReference().isPlus() ? regNew.extendBack(tDelta) : regNew.extendFront(tDelta);
                        regionChanged = true;
                    }
                    if (regionChanged) {
                        DefaultAlignedReadsData dataNew = shrinkReadsData(r, i, fDelta < 0 ? Math.abs(fDelta) : 0, tDelta < 0 ? Math.abs(tDelta) : 0, genomic, fPrime == null, tPrime == null);
                        outList.add(new ImmutableReferenceGenomicRegion<>(r.getReference(), regNew, dataNew));
                    } else {
                        factory.add(data, i);
                    }
                }
                if (factory.getDistinctSequences() > 0) {
                    outList.add(new ImmutableReferenceGenomicRegion<>(r.getReference(), r.getRegion(), factory.create()));
                }
                return EI.wrap(outList);
            }));
        }

        out.fill(it);
    }

    private static void fillFactory(AlignedReadsDataFactory factory, DefaultAlignedReadsData data, int distinct) {
        factory.setMultiplicity(data.getMultiplicity(distinct));
        int[] count = new int[data.getNumConditions()];
        for (int j = 0; j < count.length; j++) {
            count[j] = data.getCount(distinct, j);
        }
        factory.setCount(count);
        factory.setGeometry(data.getGeometryBeforeOverlap(distinct), data.getGeometryOverlap(distinct), data.getGeometryAfterOverlap(distinct));
    }

    private static DefaultAlignedReadsData shrinkReadsData(ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> read, int distinct, int shrink5p, int shrink3P, Genomic genomic, boolean add5Psoftclip, boolean add3Psoftclip) {
        DefaultAlignedReadsData data = read.getData();
        AlignedReadsDataFactory factory = new AlignedReadsDataFactory(data.getNumConditions());
        factory.start();
        factory.newDistinctSequence();
        if (add5Psoftclip) {
            GenomicRegion reg;
            if (read.getReference().isPlus()) {
                reg = new ArrayGenomicRegion(read.getRegion().getStart(), read.getRegion().getStart() + shrink5p);
            } else {
                reg = new ArrayGenomicRegion(read.getRegion().getEnd() - shrink5p, read.getRegion().getEnd());
            }
            factory.addVariation(new AlignedReadsSoftclip(true, genomic.getSequence(read.getReference(), reg), false));
            factory.addVariation(new AlignedReadsSoftclip(true, StringUtils.reverse(genomic.getSequence(read.getReference().toOppositeStrand(), reg)), true));
        }
        if (add3Psoftclip) {
            GenomicRegion reg;
            if (read.getReference().isPlus()) {
                reg = new ArrayGenomicRegion(read.getRegion().getEnd() - shrink3P, read.getRegion().getEnd());
            } else {
                reg = new ArrayGenomicRegion(read.getRegion().getStart(), read.getRegion().getStart() + shrink3P);
            }
            factory.addVariation(new AlignedReadsSoftclip(false, genomic.getSequence(read.getReference(), reg), false));
            factory.addVariation(new AlignedReadsSoftclip(false, StringUtils.reverse(genomic.getSequence(read.getReference().toOppositeStrand(), reg)), true));
        }
        for (int j = 0; j < data.getVariationCount(distinct); j++) {
            AlignedReadsVariation var = data.getVariation(distinct, j);
            if (var.isSoftclip()) {
                AlignedReadsSoftclip softclip = (AlignedReadsSoftclip)var;

                CharSequence seq;
                if (read.getReference().isPlus()) {
                    if (softclip.getPosition() == 0) {
                        GenomicRegion reg = new ArrayGenomicRegion(read.getRegion().getStart(), read.getRegion().getStart() + shrink5p);
                        if (!softclip.isFromSecondRead()) {
                            seq = softclip.getReadSequence().toString() + genomic.getSequence(read.getReference(), reg);
                        } else {
                            seq = softclip.getReadSequence().toString() + SequenceUtils.getDnaComplement(genomic.getSequence(read.getReference(), reg));
                        }
                    } else {
                        GenomicRegion reg = new ArrayGenomicRegion(read.getRegion().getEnd() - shrink3P, read.getRegion().getEnd());
                        if (!softclip.isFromSecondRead()) {
                            seq = genomic.getSequence(read.getReference(), reg) + softclip.getReadSequence().toString();
                        } else {
                            seq = SequenceUtils.getDnaComplement(genomic.getSequence(read.getReference(), reg)) + softclip.getReadSequence().toString();
                        }
                    }
                } else {
                    if (softclip.getPosition() == 0) {
                        GenomicRegion reg = new ArrayGenomicRegion(read.getRegion().getEnd() - shrink5p, read.getRegion().getEnd());
                        if (!softclip.isFromSecondRead()) {
                            seq = softclip.getReadSequence().toString() + genomic.getSequence(read.getReference(), reg);
                        } else {
                            seq = softclip.getReadSequence().toString() + SequenceUtils.getDnaComplement(genomic.getSequence(read.getReference(), reg));
                        }
                    } else {
                        GenomicRegion reg = new ArrayGenomicRegion(read.getRegion().getStart(), read.getRegion().getStart() + shrink3P);
                        if (!softclip.isFromSecondRead()) {
                            seq = genomic.getSequence(read.getReference(), reg) + softclip.getReadSequence().toString();
                        } else {
                            seq = SequenceUtils.getDnaComplement(genomic.getSequence(read.getReference(), reg)) + softclip.getReadSequence().toString();
                        }
                    }
                }

                factory.addVariation(new AlignedReadsSoftclip(softclip.getPosition()==0, seq, softclip.isFromSecondRead()));
            }
            else if (var.isMismatch()) {
                AlignedReadsMismatch mismatch = (AlignedReadsMismatch)var;
                if (mismatch.getPosition() - shrink5p < 0 || mismatch.getPosition() >= read.getRegion().getTotalLength()-shrink3P) {
                    continue;
                }
                factory.addVariation(new AlignedReadsMismatch(mismatch.getPosition()-shrink5p, mismatch.getReferenceSequence(), mismatch.getReadSequence(), mismatch.isFromSecondRead()));
            }
            else if (var.isInsertion()) {
                AlignedReadsInsertion insertion = (AlignedReadsInsertion)var;
                if (insertion.getPosition() - shrink5p < 0 || insertion.getPosition() >= read.getRegion().getTotalLength()-shrink3P) {
                    continue;
                }
                factory.addVariation(new AlignedReadsInsertion(insertion.getPosition()-shrink5p, insertion.getReadSequence(), insertion.isFromSecondRead()));
            }
            else if (var.isDeletion()) {
                AlignedReadsDeletion deletion = (AlignedReadsDeletion)var;
                if (deletion.getPosition() - shrink5p < 0 || deletion.getPosition() >= read.getRegion().getTotalLength()-shrink3P) {
                    continue;
                }
                factory.addVariation(new AlignedReadsDeletion(deletion.getPosition()-shrink5p, deletion.getReferenceSequence(), deletion.isFromSecondRead()));
            }
        }
        fillFactory(factory, data, distinct);
        return factory.create();
    }

    private static AlignedReadsSoftclip getSoftclip(AlignedReadsData data, int distinct, boolean threePrime) {
        for (int j = 0; j < data.getVariationCount(distinct); j++) {
            AlignedReadsVariation var = data.getVariation(distinct, j);
            if (var.isSoftclip()) {
                AlignedReadsSoftclip softclip = (AlignedReadsSoftclip)var;
                if (!threePrime && softclip.getPosition() == 0 && !softclip.isFromSecondRead()) {
                    return softclip;
                }
                if (threePrime && softclip.getPosition() > 0 && softclip.isFromSecondRead()) {
                    return softclip;
                }
            }
        }
        return null;
    }

    private static int softClipNum(DefaultAlignedReadsData data, int distinct) {
        int softClipNum = 0;
        for (int i = 0; i < data.getVariationCount(distinct); i++) {
            if (data.getVariation(distinct, i).isSoftclip()) {
                softClipNum++;
            }
        }
        return softClipNum;
    }
}
