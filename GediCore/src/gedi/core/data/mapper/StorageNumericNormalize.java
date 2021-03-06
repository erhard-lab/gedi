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
package gedi.core.data.mapper;

import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.tree.redblacktree.IntervalTree;

import java.util.function.UnaryOperator;

@GenomicRegionDataMapping(fromType=IntervalTree.class,toType=IntervalTree.class)
public class StorageNumericNormalize extends StorageNumericCompute {

	

	public StorageNumericNormalize(double[] totals) {
		super(new RpmNormalize(totals));
	}

	private static class RpmNormalize implements UnaryOperator<NumericArray> {

		private double[] totals;

		public RpmNormalize(double[] totals) {
			this.totals = totals;
		}

		@Override
		public NumericArray apply(NumericArray t) {
			NumericArray re = t.copy();
			for (int i=0; i<re.length(); i++)
				re.mult(i,1E6/totals[i]);
			return re;
		}
		
	}


	
}
