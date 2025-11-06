package gedi.iTiSS.data.dependencyTree;

import gedi.util.StringUtils;
import gedi.util.functions.EI;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class DependencyTree<T> {
    private Set<DependencyNode<T>> root;
    private T data;

    public DependencyTree() {

    }

    public DependencyTree(String depString) {
        buildTree(depString);
    }

    public void buildTree(String depString) {
        // depString example: 0:1>6x0,2:3>7x0,4:5>8x0,6:7:8>9x2
        String[] entity = StringUtils.split(depString, ',');
        Set<Integer> outNodes = new HashSet<>();
        ArrayList<DependencyNodeMeta> depNodes = EI.wrap(entity).map(e -> {
            String[] inOut = StringUtils.split(e, '>');
            String[] outS = StringUtils.split(inOut[1], 'x');
            int out = Integer.parseInt(outS[0]);
            outNodes.add(out);
            int tsrThresh = Integer.parseInt(outS[1]);
            int tissThresh = outS.length == 3 ? Integer.parseInt(outS[2]) : 1;
            int[] in = EI.wrap(StringUtils.split(inOut[0], ':')).map(Integer::parseInt).toIntArray();
            return new DependencyNodeMeta(in, out, tsrThresh, tissThresh);
        }).list();
        // fill roots
        root = new HashSet<>();
        for (int i = 0; i < depNodes.size(); i++) {
            for (int j = 0; j < depNodes.get(i).getInput().length; j++) {
                if (!outNodes.contains(depNodes.get(i).getInput()[j])) {
                    if (getNodeWithId(depNodes.get(i).getInput()[j]) != null) {
                        continue;
                    }
                    root.add(new DependencyNode<T>(depNodes.get(i).getInput()[j], null, 0, 1));
                }
            }
        }

        while (depNodes.size() > 0) {
            List<DependencyNodeMeta> nodesToRemove = new ArrayList<>();
            for (int i = 0; i < depNodes.size(); i++) {
                DependencyNodeMeta newNode = depNodes.get(i);
                List<DependencyNode<T>> inputNodes = new ArrayList<>();
                for (int j = 0; j < newNode.getInput().length; j++) {
                    DependencyNode<T> n = getNodeWithId(newNode.getInput()[j]);
                    if (n == null) {
                        break;
                    }
                    inputNodes.add(n);
                }
                if (inputNodes.size() == newNode.getInput().length) {
                    nodesToRemove.add(newNode);
                    DependencyNode<T> dependencyNode = new DependencyNode<T>(newNode.getOutput(), inputNodes, newNode.getTsrThresh(), newNode.getTissThresh());
                    for (DependencyNode<T> n : inputNodes) {
                        n.addOutput(dependencyNode);
                    }
                }
            }
            depNodes.removeAll(nodesToRemove);
        }
    }

    public DependencyNode<T> getNodeWithId(int id) {
        for (DependencyNode<T> node : root) {
            DependencyNode<T> ret = node.getNodeWithId(id);
            if (ret != null) {
                return ret;
            }
        }
        return null;
    }

    public Set<DependencyNode<T>> getLeaves() {
        Set<DependencyNode<T>> leaves = new HashSet<>();
        for (DependencyNode<T> node : root) {
            getLeavesHelper(node, leaves);
        }
        return leaves;
    }

    private void getLeavesHelper(DependencyNode<T> node, Set<DependencyNode<T>> leaves) {
        if (node.getOutNodes() == null) {
            leaves.add(node);
        } else {
            for (DependencyNode<T> nextNode : node.getOutNodes()) {
                getLeavesHelper(nextNode, leaves);
            }
        }
    }

    public Set<DependencyNode<T>> getRoot() {
        return root;
    }

    public void setData(T data) {
        this.data = data;
    }

    public T getData() {
        return data;
    }
}
