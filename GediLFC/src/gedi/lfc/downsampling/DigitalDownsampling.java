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
package gedi.lfc.downsampling;

import gedi.lfc.Downsampling;

public class DigitalDownsampling extends Downsampling {

	@Override
	protected void downsample(double[] counts) {
		for (int i=0; i<counts.length; i++)
			counts[i] = Math.min(1, counts[i]);
	}
}


