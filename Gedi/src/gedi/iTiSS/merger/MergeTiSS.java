package gedi.iTiSS.merger;

import gedi.core.reference.ReferenceSequence;
import gedi.iTiSS.TiSSMergerParameterSet;
import gedi.iTiSS.data.dependencyTree.DependencyNode;
import gedi.iTiSS.data.dependencyTree.DependencyTree;
import gedi.iTiSS.utils.TiSSUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;

import java.io.IOException;
import java.util.*;

public class MergeTiSS extends GediProgram {

    public MergeTiSS(TiSSMergerParameterSet params) {
        addInput(params.inputFiles);
        addInput(params.dependencies);
        addInput(params.prefix);
        addInput(params.pacbioIndices);
        addInput(params.minionIndices);
        addInput(params.pacbioDelta);
        addInput(params.minionDelta);
        addInput(params.priorityDatasets);
        addInput(params.gapMerge);

        addOutput(params.outFile);
    }

    @Override
    public String execute(GediProgramContext context) throws Exception {
        List<String> inputFiles = getParameters(0);
        String dependencies = getParameter(1);
        String prefix = getParameter(2);
        String pacbio = getParameter(3);
        String minion = getParameter(4);
        int pacbioDelta = getParameter(5);
        int minionDelta = getParameter(6);
        String priotiyDatasets = getParameter(7);
        int gap = getParameter(8);

        DependencyTree<Map<ReferenceSequence, List<Tsr>>> tree = mergeAllTiss(inputFiles, dependencies, prefix,
                pacbio, minion, pacbioDelta, minionDelta,priotiyDatasets, context, false, gap);

        int fileCount = 1;
        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> leaf : tree.getLeaves()) {
            System.err.println("Leave node with ID: " + leaf.getId());
            LineWriter writer = new LineOrientedFile(prefix + fileCount + "TiSS.tsv").write();
            LineWriter writer2 = new LineOrientedFile(prefix + fileCount + ".TSRs.tsv").write();
            writer.write("Reference\tTiSS\tTSR");
            writer2.write("Reference\tMaxTiSS\tTSR");
            for (int t = 0; t < tree.getRoot().size(); t++) {
                writer.write("\t" + inputFiles.get(t));
                writer2.write("\t" + inputFiles.get(t));
            }
            writer.writeLine("\tTiSS-Score\tTSR-Score");
            writer2.writeLine("\tMaxTiSS-Score\tTSR-Score");
            TiSSUtils.mergeTsrs(leaf.getData());
            for (ReferenceSequence ref : leaf.getData().keySet()) {
                for (Tsr i : leaf.getData().get(ref)) {
                    Tiss tiss = i.getSingleMaxTiSS();
                    writer2.write(ref.toPlusMinusString() + "\t" + tiss.getTissPos());
                    writer2.write("\t" + i.getStart() + "-" + i.getEnd());
                    for (int t = 0; t < tree.getRoot().size(); t++) {
                        writer2.write("\t");
                        if (i.getCalledBy().contains(t)) {
                            writer2.write("1");
                        } else {
                            writer2.write("0");
                        }
                    }
                    writer2.writeLine("\t" + tiss.getCalledBy().size() + "\t" + i.getCalledBy().size());
                    for (Tiss maxTss : i.getMaxTiss()) {
                        writer.write(ref.toPlusMinusString() + "\t" + maxTss.getTissPos());
                        writer.write("\t" + i.getStart() + "-" + i.getEnd());
                        int score = 0;
                        for (int t = 0; t < tree.getRoot().size(); t++) {
                            writer.write("\t");
                            if (maxTss.getCalledBy().contains(t)) {
                                writer.write("1");
                                score++;
                            } else{
                                writer.write("0");
                            }
                        }
                        writer.writeLine("\t" + score + "\t" + i.getCalledBy().size());
                    }
                }
            }
            writer.close();
            writer2.close();
        }
        return null;
    }

    public DependencyTree<Map<ReferenceSequence, List<Tsr>>> mergeAllTiss(List<String> inputFiles, String dependencies, String prefix, String pacbio,
                                                                          String minion, int pacbioDelta, int minionDelta, String priotiyDatasets,
                                                                          GediProgramContext context, boolean isTest, int gap) throws IOException {
        Set<Integer> pacbioIndices = extractIndices(pacbio);
        Set<Integer> minionIndices = extractIndices(minion);
        Set<Integer> prioritiesIndices = extractIndices(priotiyDatasets);

        DependencyTree<Map<ReferenceSequence, List<Tsr>>> tree = new DependencyTree<>(dependencies);
        Set<DependencyNode<Map<ReferenceSequence, List<Tsr>>>> root = tree.getRoot();

        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> node : root) {
            if (!isTest) {
                context.getLog().info("File: " + inputFiles.get(node.getId()) + ", ID: " + node.getId());
            }

            Map<ReferenceSequence, List<Tsr>> tsrs = TiSSUtils.extractTsrsFromFile(inputFiles.get(node.getId()), 1, gap, node.getId());
            for (ReferenceSequence ref : tsrs.keySet()) {
                EI.wrap(tsrs.get(ref)).forEachRemaining(t -> {
                    if (pacbioIndices.size() > 0 && pacbioIndices.contains(node.getId())) {
                        t.extendBack(pacbioDelta);
                        t.extendFron(pacbioDelta);
                    }
                    else if (minionIndices.size() > 0 && minionIndices.contains(node.getId())) {
                        if (ref.isPlus()) {
                            t.extendFron(minionDelta);
                        } else {
                            t.extendBack(minionDelta);
                        }
                    }
                });
            }
            node.setData(tsrs);
            node.setChecked(true);
        }

        Set<DependencyNode<Map<ReferenceSequence, List<Tsr>>>> leaves = tree.getLeaves();

        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> mapDependencyNode : root) {
            workNodes(mapDependencyNode, gap);
        }

        if (!allLeavesChecked(leaves)) {
            throw new IllegalStateException("All leaves should be checked until here!");
        }

        // TODO: Abstract!
        if (prioritiesIndices.size() > 0) {
            for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> leaf : leaves) {
                for (ReferenceSequence ref : leaf.getData().keySet()) {
                    List<Tsr> lst = leaf.getData().get(ref);
                    if (lst.size() <= 1) {
                        continue;
                    }
                    lst.sort(Tsr::compare);
                    List<Tsr> successives = new ArrayList<>();
                    Tsr lastChecked = lst.get(0);
                    List<Tsr> toDelete = new ArrayList<>();
                    for (int i = 1; i < lst.size(); i++) {
                        Tsr currentCheck = lst.get(i);
                        if (Math.abs(currentCheck.getStart() - lastChecked.getEnd()) == 0) {
                            if (successives.size() == 0) {
                                successives.add(lastChecked);
                            }
                            successives.add(currentCheck);
                        } else {
                            List<Tsr> keepers = new ArrayList<>();
                            for (Tsr succ : successives) {
                                for (int prio : prioritiesIndices) {
                                    if (tree.getNodeWithId(prio).getData().get(ref).contains(succ)) {
                                        keepers.add(succ);
                                    }
                                }
                            }
                            if (keepers.size() != 0) {
                                successives.removeAll(keepers);
                                toDelete.addAll(successives);
                            }
                            successives.clear();
                        }
                        lastChecked = currentCheck;
                    }
                    lst.removeAll(toDelete);
                }
            }
        }

        return tree;
    }

    private Set<Integer> extractIndices(String indexString) {
        Set<Integer> out = new HashSet<>();

        int[][] indicesTmp = TiSSUtils.extractReplicatesFromString(indexString, '_');

        if (indicesTmp.length == 1) {
            out = EI.wrap(indicesTmp[0]).set();
        } else if (indicesTmp.length > 1){
            throw new IllegalArgumentException("pacbio/minion string is only allowed to contain one character and the skip-character");
        }

        return out;
    }

    private void workNodes(DependencyNode<Map<ReferenceSequence, List<Tsr>>> node, int gap) {
        if (node.isChecked()) {
            if (node.getOutNodes() == null) {
                return;
            }
            for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> outNode : node.getOutNodes()) {
                workNodes(outNode, gap);
            }
            return;
        }
        if (!allNodesChecked(node.getInNodes())) {
            return;
        }
        int tsrThresh = node.getTsrThresh() <= 0 ? node.getInNodes().size() : node.getTsrThresh();
        int tissThresh = node.getTissThresh() <= 0 ? node.getInNodes().size() : node.getTissThresh();
        Map<ReferenceSequence, List<Tsr>> newData = new HashMap<>();
        List<Map<ReferenceSequence, List<Tsr>>> oldDatas = new ArrayList<>();
        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> inNode : node.getInNodes()) {
            oldDatas.add(inNode.getData());
        }
        for (Map<ReferenceSequence, List<Tsr>> oldData : oldDatas) {
            for (ReferenceSequence ref : oldData.keySet()) {
                newData.computeIfAbsent(ref, r -> new ArrayList<>());
                for (Tsr i : oldData.get(ref)) {
                    addDataToTsr(i, ref, oldDatas);
                    if (tissThresh != 1) {
                        i.removeTiss(tissThresh);
                        List<Tsr> breakUp = i.breakUp(gap);
                        breakUp.forEach(t -> {
                            if (t.getCalledBy().size() >= tsrThresh && !newData.get(ref).contains(t)) {
                                newData.get(ref).add(t);
                            }
                        });
                    } else if (i.getCalledBy().size() >= tsrThresh && !newData.get(ref).contains(i)) {
                        newData.get(ref).add(i);
                    }
                }
            }
        }
        TiSSUtils.mergeTsrs(newData);
        node.setChecked(true);
        node.setData(newData);
        if (node.getOutNodes() == null) {
            return;
        }
        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> outNode : node.getOutNodes()) {
            workNodes(outNode, gap);
        }
    }

    private void addDataToTsr(Tsr value, ReferenceSequence lookUpRef, List<Map<ReferenceSequence, List<Tsr>>> data) {
        Tiss ts = value.getSingleMaxTiSS();
        for (Map<ReferenceSequence, List<Tsr>> tsrMap : data) {
            List<Tsr> tsrs = tsrMap.get(lookUpRef);
            if (tsrs == null) {
                continue;
            }
            // TODO expensive for big data (n-squared)
            for (Tsr tsr : tsrs) {
                if (value.containsMaxTiSS(tsr)) {
                    Set<Tiss> containedTiss = value.getContainingMaxTiSS(tsr);
                    value.addAll(containedTiss);
                }
            }
        }
    }

    private boolean allLeavesChecked(Set<DependencyNode<Map<ReferenceSequence, List<Tsr>>>> leaves) {
        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> node : leaves) {
            if (!node.isChecked()) {
                return false;
            }
        }
        return true;
    }

    private boolean allNodesChecked(List<DependencyNode<Map<ReferenceSequence, List<Tsr>>>> nodes) {
        for (DependencyNode<Map<ReferenceSequence, List<Tsr>>> node : nodes) {
            if (!node.isChecked()) {
                return false;
            }
        }
        return true;
    }
}
