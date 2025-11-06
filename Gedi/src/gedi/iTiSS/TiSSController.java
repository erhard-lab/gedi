package gedi.iTiSS;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.iTiSS.analyzer.AnalyzeCustom;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.data.DataWrapper;
import gedi.iTiSS.modules.*;
import gedi.iTiSS.utils.AnalysisModuleType;
import gedi.iTiSS.utils.ReadType;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.util.StringUtils;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TiSSController extends GediProgram {
    private final char SKIP_CHAR = '_';

    public TiSSController(TiSSParameterSet params) {
        addInput(params.reads);
        addInput(params.genomic);
        addInput(params.wSize);
        addInput(params.iqr);
        addInput(params.minReadDens);
        addInput(params.replicates);
        addInput(params.pseudoCount);
        addInput(params.peakFCThreshold);
        addInput(params.timecourses);
        addInput(params.strandness);
        addInput(params.useAutoparam);
        addInput(params.plotMM);
        addInput(params.analyzeModuleType);
        addInput(params.pValThresh);
        addInput(params.cleanupThresh);
        addInput(params.testChromosomes);
        addInput(params.readType);
        addInput(params.useUpAndDownstream);
        addInput(params.minReadNum);

        addInput(params.prefix);

        addOutput(params.outIQR);
        addOutput(params.outTA);
        addOutput(params.outKA);
        addOutput(params.outXRN1);
        addOutput(params.outEquTrans);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        List<GenomicRegionStorage<AlignedReadsData>> reads = getParameters(0);
        Genomic genomic = getParameter(1);
        int windowSize = getParameter(2);
        double iqr = getParameter(3);
        double minReadDens = getParameter(4);
        String replicates = getParameter(5);
        double pseudoCount = getParameter(6);
        double peakFCThreshold = getParameter(7);
        String timecourses = getParameter(8);
        Strandness strandness = getParameter(9);
        boolean useMM = getParameter(10);
        boolean plotMM = getParameter(11);
        List<AnalysisModuleType> modTypes = getParameters(12);
        double pValThresh = getParameter(13);
        int cleanupThresh = getParameter(14);
        String testChromosomes = getParameter(15);
        ReadType readType = getParameter(16);
        boolean useUpAndDownstream = getParameter(17);
        int minReadNum = getParameter(18);

        String prefix = getParameter(19);

        final boolean useMultiCourse = timecourses != null && !timecourses.isEmpty();

        if (replicates == null || replicates.isEmpty()) {
            context.getLog().severe("The parameter -rep is missing!!");
            context.getLog().severe("If you are running iTiSS on a single dataset, simply add '-rep X' to your call.");
            context.getLog().severe("If you are running it on multiple datasets at once, please visit the Wiki on GitHub.");
        }

        int[][] reps = TiSSUtils.extractReplicatesFromString(replicates, SKIP_CHAR);

        Set<ReferenceSequence> testChrs = new HashSet<>();
        if (testChromosomes != null && !testChromosomes.isEmpty()) {
            for (String r : StringUtils.split(testChromosomes, ",")) {
                if (r != null && !r.isEmpty()) {
                    testChrs.add(Chromosome.obtain(r));
                }
            }
        }
        if (testChrs.size() == 0) {
            testChrs = null;
        }

        DataWrapper dataWrapper = new DataWrapper(reads, strandness, readType, testChrs);
        List<Data> data = new ArrayList<>();
        List<Data> singleLanes = new ArrayList<>(reps.length);
        List<Data> multiLanes = new ArrayList<>(reps.length);
        for (int[] rep : reps) {
            if (!useMultiCourse) {
                Data d = new Data(rep, false);
                data.add(d);
                singleLanes.add(d);
                context.getLog().info("Using non-multi data handling.");
            }
            else {
                Data d = new Data(rep, true);
                data.add(d);
                multiLanes.add(d);
                context.getLog().info("Using multi data handling.");
            }
        }
        dataWrapper.initData(genomic, data);

        AnalyzeCustom analyzer = new AnalyzeCustom();

        for (AnalysisModuleType modType : modTypes) {
            switch (modType) {
                case DENSITY:
                    context.getLog().info("Adding Density module");
                    singleLanes.forEach(d -> analyzer.addModule(new TranscriptionalActivity(pValThresh, windowSize, minReadDens, cleanupThresh, minReadNum, useMM, prefix, d, "DENSITY")));
                    break;
                case KINETIC:
                    context.getLog().info("Adding Kinetic module");
                    multiLanes.forEach(d -> analyzer.addModule(new KineticActivity(windowSize, pValThresh, pseudoCount, cleanupThresh, useMM, useUpAndDownstream, minReadNum, prefix, d, "KINETIC")));
                    break;
                case SPARSE_PEAK:
                    context.getLog().info("Adding sparse peak module");
                    singleLanes.forEach(d -> analyzer.addModule(new DRnaModule(windowSize, pseudoCount, peakFCThreshold, cleanupThresh, minReadNum, d, useMM, "SPARSE_PEAK")));
                    break;
                case DENSE_PEAK:
                    context.getLog().info("Adding dense peak module");
                    singleLanes.forEach(d -> analyzer.addModule(new CRnaModule(iqr, pseudoCount, windowSize, cleanupThresh, minReadNum, d, useMM, "DENSE_PEAK")));
                    break;
                case EQUAL_TRANSCRIPTION:
                    context.getLog().info("Adding equal transcription module");
                    singleLanes.forEach(d -> analyzer.addModule(new EqualTranscriptionModule("EQUAL_TRANSCRIPTION", d)));
                    break;
                default:
                    throw new IllegalArgumentException("AnalyzeModuleType is unrecognized: " + modType);
            }
        }

        analyzer.startAnalyzing(dataWrapper, genomic);

        context.getLog().info("Analyzation modules finished");
        context.getLog().info("Writing final file(s)");

        for (AnalysisModuleType modType : modTypes) {
            switch (modType) {
                case DENSITY:
                    if (useMM) {
                        analyzer.calculateMLResults(TranscriptionalActivity.class, prefix, plotMM);
                    }
                    analyzer.writeOutResults(new LineOrientedFile(getOutputFile(1).getPath()), TranscriptionalActivity.class);
                    break;
                case KINETIC:
                    if (useMM) {
                        analyzer.calculateMLResults(KineticActivity.class, prefix, plotMM);
                    }
                    analyzer.writeOutResults(new LineOrientedFile(getOutputFile(2).getPath()), KineticActivity.class);
                    break;
                case SPARSE_PEAK:
                    if (useMM) {
                        analyzer.calculateMLResults(DRnaModule.class, prefix, plotMM);
                    }
                    analyzer.writeOutResults(new LineOrientedFile(getOutputFile(3).getPath()), DRnaModule.class);
                    break;
                case DENSE_PEAK:
                    if (useMM) {
                        analyzer.calculateMLResults(CRnaModule.class, prefix, plotMM);
                    }
                    analyzer.writeOutResults(new LineOrientedFile(getOutputFile(0).getPath()), CRnaModule.class);
                    break;
                case EQUAL_TRANSCRIPTION:
                    analyzer.writeOutResults(new LineOrientedFile(getOutputFile(4).getPath()), EqualTranscriptionModule.class);
                    break;
                default:
                    throw new IllegalArgumentException("AnalysisModuleType is unrecognized: " + modType);
            }
        }

        return null;
    }
}
