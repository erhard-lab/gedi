package gedi.iTiSS.utils.nonGediCompatible;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.region.GenomicRegionStorage;
import gedi.iTiSS.TiSSParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.util.List;

public class ConvertBam extends GediProgram {
    public ConvertBam(TiSSParameterSet params) {
        addInput(params.bamFiles);
        addInput(params.prefix);

        addOutput(params.reads);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        List<String> bamFiles = getParameters(0);
        String prefix = getParameter(1);

        if (bamFiles.size() == 0) {
            context.getLog().info("No BAM files provided. Using CIT-file.");
            return null;
        }
        context.getLog().info("Converting BAMs into CIT-files. This might take a while depending on the BAM size");
        context.getLog().info("Keep in mind that they need to be sorted and indexed using samtools prior to iTiSS!");
        String readsOutputFile = prefix + "BAM2CIT_convertedReads.cit";

        GenomicRegionStorage<DefaultAlignedReadsData> cit = BamConverter.convertBam(bamFiles, readsOutputFile);

        setOutput(0, cit);
        return null;
    }
}
