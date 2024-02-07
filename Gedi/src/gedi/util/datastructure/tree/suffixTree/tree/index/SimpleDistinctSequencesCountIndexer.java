package gedi.util.datastructure.tree.suffixTree.tree.index;

import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.suffixTree.construction.UkkonenSuffixTreeBuilder;
import gedi.util.datastructure.tree.suffixTree.tree.GeneralizedSuffixTreeIndex;
import gedi.util.datastructure.tree.suffixTree.tree.SuffixTree;
import gedi.util.datastructure.tree.suffixTree.tree.traversal.DfsDownAndUpTraverser;
import gedi.util.datastructure.tree.suffixTree.tree.traversal.Traverser;

public class SimpleDistinctSequencesCountIndexer extends AbstractIntIndexer {

	private GeneralizedSuffixTreeIndex ukk;
	
	
	public SimpleDistinctSequencesCountIndexer(GeneralizedSuffixTreeIndex ukk) {
		this.ukk = ukk;
	}

	@Override
	protected void createNew_internal(SuffixTree tree, int[] index) {
		IntArrayList[] sets = new IntArrayList[index.length];
		
		DfsDownAndUpTraverser t = new DfsDownAndUpTraverser(tree,tree.getRoot().getNode());
		while (t.hasNext()) {
			int node = t.nextInt();
			
			int direction = t.getDirection();
			if (direction==Traverser.UP) {
				if (sets[node]==null) sets[node] = new IntArrayList(2);
				if (sets[t.getPrevious()]==null) {
					int pos = tree.getIntAttributes(UkkonenSuffixTreeBuilder.SUFFIX_LINK_NAME)[t.getPrevious()];
					int si = ukk.getSequenceIndexForGeneralizedPosition(pos);
					index[t.getPrevious()] = 1; // to store sequence index: si;
					sets[node].add(si);
				} else {
					sets[t.getPrevious()].sort();
					sets[t.getPrevious()].unique();
					index[t.getPrevious()] = sets[t.getPrevious()].size();
					sets[node].addAll(sets[t.getPrevious()]);
					sets[t.getPrevious()] = null;
				}
			}
		}
		sets[0].sort();
		sets[0].unique();
		index[0] = sets[0].size();
	}

	@Override
	public String name() {
		return "DistinctSequencesCount";
	}

}