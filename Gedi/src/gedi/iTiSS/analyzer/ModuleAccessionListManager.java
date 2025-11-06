package gedi.iTiSS.analyzer;

import gedi.iTiSS.modules.ModuleBase;
import gedi.iTiSS.utils.scheduling.AccessionList;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class ModuleAccessionListManager {
    private Map<ModuleBase, AccessionList> moduleAccessList;

    public ModuleAccessionListManager(Set<ModuleBase> modules, int accessCount) {
        this.moduleAccessList = new HashMap<>();
        for (ModuleBase module : modules) {
            moduleAccessList.put(module, new AccessionList(accessCount));
        }
    }

    public synchronized int accessNextFree(ModuleBase module) {
        int next = moduleAccessList.get(module).getNextNeededReset();
        if (next == -1) {
            throw new IllegalStateException("Nothing to access. Are all the modules already finished?");
        }
//        int looped = next;
//        while (isUsed(module, next)) {
//            next = moduleAccessList.get(module).getNextNeeded();
//            if (next == looped) {
//                return -1;
//            }
//        }
        while (isUsed(module, next)) {
            next = moduleAccessList.get(module).getNextNeeded();
            if (next == -1) {
                return next;
            }
        }
        moduleAccessList.get(module).startAccess(next);
        return next;
    }

    public synchronized void finishAccess(ModuleBase module, int index) {
        moduleAccessList.get(module).stopAccess(index);
    }

    public boolean allFinished(ModuleBase module) {
        return moduleAccessList.get(module).allFinished();
    }

    private boolean isUsed(ModuleBase module, int access) {
        for (ModuleBase moduleBase : moduleAccessList.keySet()) {
            if (moduleBase.hashCode() == module.hashCode()) {
                continue;
            }
            if (!moduleAccessList.get(moduleBase).isFree(access)) {
                return true;
            }
        }
        return false;
    }
}
