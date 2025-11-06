package gedi.iTiSS.utils.nonGediCompatible;

import gedi.core.genomic.Genomic;
import gedi.iTiSS.TiSSParameterSet;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.util.ArrayList;
import java.util.List;

public class GenomicCreate extends GediProgram {

    public GenomicCreate(TiSSParameterSet params) {
//        addInput(params.genomicName);
//        addInput(params.genomicLength);
        addInput(params.chromSizes);
        addInput(params.bamFiles);
        addInput(params.prefix);

        addOutput(params.genomic);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
//        List<String> genomicName = getParameters(0);
//        List<Integer> genomicLength = getParameters(1);
        String chromSizesFile = getParameter(0);

        if (chromSizesFile == null || chromSizesFile.isEmpty()) {
            context.getLog().info("No genomic length and name provided. Using indexed genomic object.");
            return null;
        }

        List<String> genomicName = new ArrayList<>();
        List<Integer> genomicLength = new ArrayList<>();
        EI.lines(chromSizesFile).forEachRemaining(l -> {
            String[] split = StringUtils.split(l, "\t");
            genomicName.add(split[0]);
            genomicLength.add(Integer.parseInt(split[1]));
        });

        Genomic genomic = SyntheticGenome.createSyntheticGenomeOnlyTranscripts(genomicName, genomicLength);

//        genomic.iterateReferenceSequences().forEachRemaining(r -> context.getLog().info(r.toPlusMinusString() + "\t" + genomic.getLength(r.toPlusMinusString())));

        setOutput(0, genomic);
        return null;
    }
}

