/**
 * 
 *    Copyright 2017 Florian Erhard
 *
 *   Licensed under the Apache License, Version 2.0 (the "License");
 *   you may not use this file except in compliance with the License.
 *   You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *   See the License for the specific language governing permissions and
 *   limitations under the License.
 * 
 */
package gedi.util.datastructure.tree.suffixTree.tree.index;

import gedi.util.datastructure.tree.suffixTree.tree.SuffixTree;
import gedi.util.datastructure.tree.suffixTree.tree.traversal.DfsDownAndUpTraverser;
import gedi.util.datastructure.tree.suffixTree.tree.traversal.Traverser;

public class LeavesCountIndexer extends AbstractIntIndexer {

	@Override
	protected void createNew_internal(SuffixTree tree, int[] index) {
		DfsDownAndUpTraverser t = new DfsDownAndUpTraverser(tree,tree.getRoot().getNode());
		int prev = -2;
		while (t.hasNext()) {
			int node = t.nextInt();
			int direction = t.getDirection();
			if (node==prev && direction==Traverser.UP)
				index[t.getPrevious()]=1;
			if (node>=0 && direction==Traverser.UP)
				index[node] += index[t.getPrevious()];
			if (direction==Traverser.DOWN)
				prev = t.getPrevious();
		}
	}

	@Override
	public String name() {
		return "IndexLeavesCount";
	}

}
