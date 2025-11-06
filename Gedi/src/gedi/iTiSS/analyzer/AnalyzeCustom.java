package gedi.iTiSS.analyzer;

import gedi.iTiSS.modules.ModuleBase;

public class AnalyzeCustom extends AnalyzerBase{
    public AnalyzeCustom() {

    }

    public void addModule(ModuleBase module) {
        modules.add(module);
    }
}
