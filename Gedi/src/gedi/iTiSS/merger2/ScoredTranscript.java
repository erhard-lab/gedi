package gedi.iTiSS.merger2;

import executables.GenomicUtils;
import gedi.core.data.annotation.Transcript;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.iTiSS.utils.GenomicRegionUtils;

public class ScoredTranscript {
    private int totalDelta;
    private ImmutableReferenceGenomicRegion<Transcript> transcript;
    private TsrFileEntry tss;
    private TsrFileEntry tts;

    public ScoredTranscript(int totalDelta, ImmutableReferenceGenomicRegion<Transcript> transcript, TsrFileEntry tss, TsrFileEntry tts) {
        this.totalDelta = totalDelta;
        this.transcript = transcript;
        this.tss = tss;
        this.tts = tts;
    }

    public int getTss() {
        return tss == null ? GenomicRegionPosition.FivePrime.position(transcript) : tss.getMaxReadCountTissPos();
    }

    public int getTts() {
        return tts == null ? GenomicRegionPosition.ThreePrime.position(transcript) : tts.getMaxReadCountTissPos();
    }

    public ImmutableReferenceGenomicRegion<Transcript> getTranscript() {
        return transcript;
    }

    public boolean adjustTranscriptRegion() {
        GenomicRegion newTranscriptRegion = cutOutPotentialIntrons(transcript.getRegion(), getTss()); // these two are only executed
        newTranscriptRegion = cutOutPotentialIntrons(newTranscriptRegion, getTts()); // when allowIntrons is set to true
        int left = transcript.getReference().isPlus() ? getTss() : getTts();
        int right = transcript.getReference().isPlus() ? getTts() : getTss();
        if (left >= newTranscriptRegion.getEnd() || right < newTranscriptRegion.getStart() || right <= left) {
            return false;
        }
        newTranscriptRegion = newTranscriptRegion.extendFront(-newTranscriptRegion.induceMaybeOutside(left));
        newTranscriptRegion = newTranscriptRegion.extendBack(newTranscriptRegion.induceMaybeOutside(right)-newTranscriptRegion.getTotalLength());
        transcript = transcript.toMutable().setRegion(newTranscriptRegion).toImmutable();
        return true;
    }

    private static GenomicRegion cutOutPotentialIntrons(GenomicRegion region, int pos) {
        int intron = GenomicRegionUtils.getIntronIndex(region, pos);

        if (intron == -1) {
            return region;
        }

        return GenomicRegionUtils.removeIntron(region, intron);
    }

    public int getTotalDelta() {
        return totalDelta;
    }

    public GenomicRegion getRegionForHashing() {
        if (transcript.getReference().isPlus()) {
            if (tss == null) {
                return new ArrayGenomicRegion(transcript.getRegion().getEnd()-1, transcript.getRegion().getEnd());
            } else if (tts == null) {
                return new ArrayGenomicRegion(transcript.getRegion().getStart(), transcript.getRegion().getStart()+1);
            } else {
                return new ArrayGenomicRegion(transcript.getRegion().getStart(), transcript.getRegion().getEnd());
            }
        } else {
            if (tss == null) {
                return new ArrayGenomicRegion(transcript.getRegion().getStart(), transcript.getRegion().getStart()+1);
            } else if (tts == null) {
                return new ArrayGenomicRegion(transcript.getRegion().getEnd()-1, transcript.getRegion().getEnd());
            } else {
                return new ArrayGenomicRegion(transcript.getRegion().getStart(), transcript.getRegion().getEnd());
            }
        }
    }
}
