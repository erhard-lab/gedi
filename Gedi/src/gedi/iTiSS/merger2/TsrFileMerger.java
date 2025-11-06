package gedi.iTiSS.merger2;

import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class TsrFileMerger extends GediProgram {
    public TsrFileMerger(TissMerger2ParameterSet params) {
        addInput(params.prefix);
        addInput(params.inTissFiles);
        addInput(params.gap);
        addInput(params.extension);
        addInput(params.minScore);
        addInput(params.inTsrFiles);
        addInput(params.keepFileIds);

        addOutput(params.outFile);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        String prefix = getParameter(0);
        List<String> inTissFiles = getParameters(1);
        int gap = getParameter(2);
        int ext = getParameter(3);
        int minScore = getParameter(4);
        List<String> inTsrFiles = getParameters(5);
        boolean keepFileIds = getParameter(6);

        List<TsrFile> tsrFiles = readTissFiles(inTissFiles, gap);
        tsrFiles.addAll(readTsrFiles(inTsrFiles, keepFileIds));

        TsrFile merged = merge(tsrFiles, ext);
        merged = merged.filter(minScore);

        String outPath = getOutputFile(0).getPath();

        ExtendedIterator<String> originNamesIt = EI.wrap(convertFilePathsToFileNames(inTissFiles));
        if (keepFileIds) {
            originNamesIt = originNamesIt.chain(EI.wrap(getTsrFileNameAssociations(inTsrFiles)));
        } else {
            originNamesIt = originNamesIt.chain(EI.wrap(convertFilePathsToFileNames(inTsrFiles)));
        }

        String[] originNames = originNamesIt.toArray(new String[] {});
        merged.writeToFile(outPath, originNames);

        return null;
    }

    public static List<String> convertFilePathsToFileNames(List<String> paths) {
        return EI.wrap(paths).map(f -> Paths.get(f).getFileName().toString()).list();
    }

    public static List<String> getTsrFileNameAssociations(List<String> tsrFiles) throws IOException {
        List<String> fileNames = new ArrayList<>();
        for (String tsrFile : tsrFiles) {
            fileNames.addAll(EI.wrap(StringUtils.split(StringUtils.splitField(EI.lines(tsrFile).next(), " ", 1), ",")).list());
        }
        return fileNames;
    }

    public List<TsrFile> readTissFiles(List<String> inTissFiles, int gap) {
        List<TsrFile> tsrFiles = new ArrayList<>();
        for (int i = 0; i < inTissFiles.size(); i++) {
            tsrFiles.add(TissFile.loadFromFile(inTissFiles.get(i)).reduceToTsr(gap));
        }
        return tsrFiles;
    }

    public List<TsrFile> readTsrFiles(List<String> inTsrFiles, boolean keepFileIds) throws IOException {
        List<TsrFile> tsrFiles = new ArrayList<>();
        for (int i = 0; i < inTsrFiles.size(); i++) {
            TsrFile tsrFile = TsrFile.loadFromFile(inTsrFiles.get(i));
            if (!keepFileIds) {
                tsrFile.resetFileAssociations(Paths.get(inTsrFiles.get(i)).getFileName().toString().hashCode());
            }
            tsrFiles.add(tsrFile);
        }
        return tsrFiles;
    }

    public TsrFile merge(List<TsrFile> tissFiles, final int extension) {
        TsrFile finalFile = tissFiles.get(0);
        for (int i = 1; i < tissFiles.size(); i++) {
            finalFile = finalFile.merge(tissFiles.get(i), extension);
        }
        return finalFile;
    }
}
