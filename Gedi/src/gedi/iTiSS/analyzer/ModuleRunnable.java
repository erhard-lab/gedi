package gedi.iTiSS.analyzer;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.data.DataWrapper;
import gedi.iTiSS.modules.ModuleBase;
import gedi.iTiSS.utils.multithreading.NotifyOnFinishedRunnable;
import gedi.util.datastructure.array.NumericArray;

public class ModuleRunnable extends NotifyOnFinishedRunnable {

    private ModuleBase module;
    private ReferenceSequence ref;
    private int access;
    private DataWrapper dataWrapper;
    private int refLength;

    public ModuleRunnable(ModuleBase module) {
        this.module = module;
    }

    public void init(DataWrapper data, int access, int refLength, ReferenceSequence ref) {
        this.dataWrapper = data;
        this.access = access;
        this.refLength = refLength;
        this.ref = ref;
    }

    @Override
    public void doRun() {
//        System.err.println("[" + Thread.currentThread().getName() + "] " + module.getModuleName() + "-module starts searching for TiSS in: " + ref.toPlusMinusString());
        NumericArray[] data;
        do {
            data = dataWrapper.startAccessingData2(module.getLane(), ref, refLength);
            if (data == null) {
                try {
//                    System.err.println(Thread.currentThread().getName() + " is sleeping");
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        } while (data == null);
        module.findTiSS(data, ref);
    }

    public ModuleBase getModule() {
        return module;
    }

    public int getAccess() {
        return access;
    }
}
