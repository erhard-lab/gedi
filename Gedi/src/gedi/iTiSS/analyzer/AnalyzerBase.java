package gedi.iTiSS.analyzer;

import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.DataWrapper;
import gedi.iTiSS.modules.KineticActivity;
import gedi.iTiSS.modules.ModuleBase;
import gedi.util.dynamic.impl.DoubleDynamicObject;
import gedi.util.dynamic.impl.IntDynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutablePair;
import gedi.util.r.RRunner;

import java.io.IOException;
import java.util.*;

public abstract class AnalyzerBase{
    protected List<ModuleBase> modules = new ArrayList<>();

    public final void startAnalyzing(DataWrapper dataWrapper, Genomic genomic) {
        //TODO: Make nice
        Map<Integer, Set<ModuleBase>> pooledModules = modulePooling();
        List<ModuleScheduler> schedulers = new ArrayList<>();
        for (Integer key : pooledModules.keySet()) {
            schedulers.add(new ModuleScheduler(pooledModules.get(key), genomic, dataWrapper));
        }
        System.err.println("Number of schedulers: " + schedulers.size());
        List<Thread> threads = new ArrayList<>();
        for (ModuleScheduler scheduler : schedulers) {
            threads.add(new Thread(scheduler::run));
        }
        for(Thread thread : threads) {
            thread.start();
        }
        for (Thread thread : threads) {
            try {
                thread.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    private Map<Integer, Set<ModuleBase>> modulePooling() {
        Map<Integer, Set<ModuleBase>> pooledModules = new HashMap<>();
        for (ModuleBase module : modules) {
            int hash = module.laneHashCode();
            if (!pooledModules.containsKey(hash)) {
                pooledModules.put(hash, new HashSet<>());
            }
            pooledModules.get(hash).add(module);
        }
        return pooledModules;
    }

    public void calculateMLResults(Class<? extends ModuleBase> type, String prefix, boolean plot) throws IOException {
        List<ModuleBase> modules = EI.wrap(this.modules).filter(type::isInstance).list();

        for (ModuleBase module : modules) {
            module.calculateMLResults(prefix, plot);
        }
    }

    public int writeOutResults(LineOrientedFile file, Class<? extends ModuleBase> type) throws IOException {
        LineWriter writer = file.write();

        List<ModuleBase> modules = EI.wrap(this.modules).filter(type::isInstance).list();

        if (modules.size() == 0) {
            return 0;
        }

        int tissCounter = 0;

        ModuleBase module = modules.get(0);
        writer.write("Reference (Strand)\tPosition");
        if (module.getResultsNew().keySet().size() > 0) {
            // TODO: This can be null if no TiSS were found (in that case however, there are other problems with the dataset...)
            ReferenceSequence ref = EI.wrap(module.getResultsNew().keySet()).filter(k -> module.getResultsNew().get(k).keySet().size() > 0).next();
            for (ModuleBase mod : modules) {
                if (mod.getResultsNew().get(ref) != null) {
                    if (mod.getResultsNew().get(ref).keySet().size() > 0) {
                        Integer pos = mod.getResultsNew().get(ref).keySet().iterator().next();
                        Set<String> headers = mod.getResultsNew().get(ref).get(pos).keySet();
                        for (String header : headers) {
                            writer.write("\t" + header);
                        }
                    }
                }
            }
        }
        writer.writeLine();
        for (ReferenceSequence ref : module.getResultsNew().keySet()) {
            Map<Integer, Map<String, Double>> tiss = module.getResultsNew().get(ref);
            for (Integer i : tiss.keySet()) {
                if (wasPeakFoundInAllLanes(i, ref, modules)) {
                    writer.write(ref.toPlusMinusString() + "\t" + i.toString());
                    for (ModuleBase mod : modules) {
                        for (String head : mod.getResultsNew().get(ref).get(i).keySet()) {
                            writer.write("\t" + mod.getResultsNew().get(ref).get(i).get(head));
                        }
                    }
                    writer.writeLine();
                    tissCounter++;
                }
            }
        }

        writer.close();

        return tissCounter;
    }

    public void writeOutMmData(LineOrientedFile file, Class<? extends ModuleBase> type) throws IOException {
        LineWriter writer = file.write();

        List<ModuleBase> modules = EI.wrap(this.modules).filter(type::isInstance).list();

        if (modules.size() == 0) {
            return;
        }

        writer.writeLine("Ref\tValue");
        for (ModuleBase module : modules) {
            for (ReferenceSequence ref : module.getAllMmData().keySet()) {
                for (double d : module.getAllMmData().get(ref)) {
                    writer.writeLine(ref.toPlusMinusString() + "\t" + d);
                }
            }
        }

        writer.close();
    }

    public void writeOutThresholds(LineOrientedFile file, Class<? extends ModuleBase> type) throws IOException {
        LineWriter writer = file.write();

        List<ModuleBase> modules = EI.wrap(this.modules).filter(type::isInstance).list();

        if (modules.size() == 0) {
            return;
        }

        writer.writeLine("Ref\tThreshold");
        for (ModuleBase module : modules) {
            for (ReferenceSequence ref : module.getThresholds().keySet()) {
                writer.writeLine(ref.toPlusMinusString() + "\t" + module.getThresholds().get(ref));
            }
        }
        writer.close();
    }

    public void plotMlData(String dataFilePath, String threshFilePath, String prefix, int colNumber) throws IOException {
        RRunner r = new RRunner(prefix+".plotMMdata.R");
        r.set("dataFile",dataFilePath);
        r.set("threshFile",threshFilePath);
        r.set("pdfFile",prefix+".thresholds.pdf");
        r.set("prefix", prefix);
        r.set("customThresh", new DoubleDynamicObject(2));
        r.set("colNumber", new IntDynamicObject(colNumber));
        r.addSource(getClass().getResourceAsStream("/resources/plotMM.R"));
        r.run(false);
    }

    private boolean wasPeakFoundInAllLanes(Integer i, ReferenceSequence ref, List<ModuleBase> modulesToCheck) {
        for (int j = 1; j < modulesToCheck.size(); j++) {
            if (!modulesToCheck.get(j).getResultsNew().get(ref).keySet().contains(i)) {
                return false;
            }
        }
        return true;
    }
}
