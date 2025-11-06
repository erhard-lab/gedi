package gedi.iTiSS.bedConverter;

import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.iTiSS.TiSStoBedParameterSet;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.tsv.formats.BedEntry;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.IOException;

public class TsrToBed extends GediProgram {

    public TsrToBed(TiSStoBedParameterSet params) {
        addInput(params.prefix);
        addInput(params.tsrFile);
        addInput(params.minScore);

        addOutput(params.bedOut);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        String prefix = getParameter(0);
        String tsrFilePath = getParameter(1);
        int minScore = getParameter(2);

        LineWriter bedWriter = new LineOrientedFile(getOutputFile(0).getPath()).write();

        EI.lines(tsrFilePath).skip(1).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            int score = Integer.parseInt(split[split.length-1]);
            if (score < minScore) return;
            ReferenceSequence ref = Chromosome.obtain(split[0]);
            GenomicRegion tsrRegion = GenomicRegion.parse(split[2]);
            int maxTissPos = Integer.parseInt(split[1]);
            GenomicRegion maxTiSSregion = new ArrayGenomicRegion(maxTissPos, maxTissPos + 1);
            ImmutableReferenceGenomicRegion<ScoreNameAnnotation> tsr = new ImmutableReferenceGenomicRegion<>(ref,
                    tsrRegion, new ScoreNameAnnotation("TSR", 1));
            ImmutableReferenceGenomicRegion<ScoreNameAnnotation> tiss = new ImmutableReferenceGenomicRegion<>(ref,
                    maxTiSSregion, new ScoreNameAnnotation("TiSS", 0));
            BedEntry tsrBe = new BedEntry(tsr);
            BedEntry tissBe = new BedEntry(tiss);
            try {
                bedWriter.writeLine(tsrBe.toString());
                bedWriter.writeLine(tissBe.toString());
            } catch (IOException e) {
                e.printStackTrace();
            }
        });

        bedWriter.close();
        return null;
    }
}
