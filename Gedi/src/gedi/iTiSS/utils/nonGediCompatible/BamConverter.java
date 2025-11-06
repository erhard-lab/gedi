package gedi.iTiSS.utils.nonGediCompatible;

import executables.Bam2CIT;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.region.GenomicRegionStorage;

import java.io.IOException;
import java.util.List;

public class BamConverter {
    public static GenomicRegionStorage<DefaultAlignedReadsData> convertBam(List<String> bamFiles, String outputPath) throws IOException {
        String[] bam2CitArgs = new String[bamFiles.size()+1];
        bam2CitArgs[0] = outputPath;
//        bam2CitArgs[0] = "";
        for (int i = 1; i < bamFiles.size()+1; i++) {
            bam2CitArgs[i] = bamFiles.get(i-1);
        }
        Bam2CIT.main(bam2CitArgs);
        GenomicRegionStorage<DefaultAlignedReadsData> cit = new CenteredDiskIntervalTreeStorage<>(outputPath, DefaultAlignedReadsData.class);
        return cit;
    }
}
