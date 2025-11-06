package gedi.iTiSS.analyzer;

import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.Data;
import gedi.iTiSS.data.DataWrapper;
import gedi.iTiSS.modules.ModuleBase;
import gedi.iTiSS.utils.multithreading.RunnableFinishedListener;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.functions.EI;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class ModuleScheduler implements RunnableFinishedListener {
//    private Data moduleLanes;
    private Map<ModuleBase, Data> moduleLanes;
    private ModuleAccessionListManager moduleAccessionListManager;
    private Map<ModuleBase, ModuleRunnable> moduleRunnables;
    private ReferenceSequence[] refs;
    private Genomic genomic;
    private DataWrapper dataWrapper;

    public ModuleScheduler(Set<ModuleBase> modules, Genomic genomic, DataWrapper dataWrapper) {
//        System.err.println("ModuleScheduler created with " + modules.size() + " modules.");
        this.dataWrapper = dataWrapper;
        this.genomic = genomic;
        this.refs = EI.wrap(dataWrapper.getLoadedChromosomes()).toArray(new ReferenceSequence[0]);
//        this.moduleLanes = EI.wrap(modules).next().getLane();
        this.moduleLanes = EI.wrap(modules).toMap(new HashMap<ModuleBase, Data>(), m -> m, ModuleBase::getLane);
        this.moduleAccessionListManager = new ModuleAccessionListManager(modules, refs.length);
        init(modules);
    }

    private void init(Set<ModuleBase> modules) {
        moduleRunnables = new HashMap<>();
        for (ModuleBase module : modules) {
            ModuleRunnable moduleRunnable = new ModuleRunnable(module);
            moduleRunnables.put(module, moduleRunnable);
            moduleRunnable.addListener(this);
        }
    }

    public void run() {
        for (ModuleBase module : moduleRunnables.keySet()) {
            run(module);
        }
        while (!allModulesFinished()) {
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }

    private boolean allModulesFinished() {
        for (ModuleBase module : moduleRunnables.keySet()) {
            if (!moduleAccessionListManager.allFinished(module)) {
                return false;
            }
        }
        return true;
    }

    private void run(ModuleBase module) {
        if (moduleAccessionListManager.allFinished(module)) {
            System.err.println(module.getModuleName() + " finished all its analysis.");
            return;
        }
        int nextAccess = moduleAccessionListManager.accessNextFree(module);
        if (nextAccess < 0) {
            ModuleWaiting moduleWaiting = new ModuleWaiting(module, 1000);
            moduleWaiting.addListener(this);
            new Thread(moduleWaiting).start();
            return;
        }
        ModuleRunnable moduleRunnable = moduleRunnables.get(module);
        moduleRunnable.init(dataWrapper, nextAccess, genomic.getLength(refs[nextAccess].toPlusMinusString()), refs[nextAccess]);
        Thread thread = new Thread(moduleRunnable);
        thread.start();
    }

    @Override
    public void notifyRunnableFinished(Runnable runnable) {
        if (runnable instanceof ModuleRunnable) {
            ModuleRunnable moduleRunnable = (ModuleRunnable) runnable;
            System.err.println("[" + Thread.currentThread().getName() + "] " + moduleRunnable.getModule().getModuleName() + " finished for: " + refs[moduleRunnable.getAccess()].toPlusMinusString());
            dataWrapper.finishAccessingData(moduleLanes.get(moduleRunnable.getModule()), refs[moduleRunnable.getAccess()]);
            System.err.println("[" + Thread.currentThread().getName() + "] " + moduleRunnable.getModule().getModuleName() + " released data for: " + refs[moduleRunnable.getAccess()].toPlusMinusString());
            moduleAccessionListManager.finishAccess(moduleRunnable.getModule(), moduleRunnable.getAccess());
            run(moduleRunnable.getModule());
        } else if (runnable instanceof ModuleWaiting) {
            ModuleWaiting moduleWaiting = (ModuleWaiting) runnable;
            System.err.print("[" + Thread.currentThread().getName() + "] " + moduleWaiting.getModule().getModuleName() + ": Waiting thread finished waiting.\r");
            ModuleBase module = moduleWaiting.getModule();
            run(module);
         } else {
            System.err.println("We should never arrive here.");
            throw new IllegalStateException("Unsafe territory. Did you forget to implement a runnable-subclass catch?");
        }
    }
}
