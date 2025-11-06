package gedi.iTiSS.merger;

import gedi.core.reference.ReferenceSequence;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Tiss2Tsr {
    private Map<ReferenceSequence, List<Tsr>> tsrsMap;

    public Tiss2Tsr(Map<ReferenceSequence, List<Tiss>> tiss, int gap) {
        tsrsMap = new HashMap<>();
        for (ReferenceSequence ref : tiss.keySet()) {
            List<Tiss> tissLst = tiss.get(ref);
            tissLst.sort(Tiss::compare);
            List<Tsr> tsrs = new ArrayList<>();
            if (tissLst.size() == 0) {
                tsrsMap.put(ref, tsrs);
            }
            Tsr tsr = new Tsr(tissLst.get(0));
            for (int i = 1; i < tissLst.size(); i++) {
                int lastTissPos = tsr.getEnd();
                Tiss currentTiss = tissLst.get(i);
                if (currentTiss.getTissPos()-lastTissPos <= gap) {
                    tsr.add(currentTiss);
                } else {
                    tsrs.add(tsr);
                    tsr = new Tsr(currentTiss);
                }
            }
            tsrs.add(tsr);
            tsrsMap.put(ref, tsrs);
        }
    }

    public Map<ReferenceSequence, List<Tsr>> getTsrsMap() {
        return tsrsMap;
    }
}
