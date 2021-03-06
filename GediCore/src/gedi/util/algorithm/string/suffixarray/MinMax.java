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
package gedi.util.algorithm.string.suffixarray;

/**
 * Holder for minimum and maximum.
 * 
 * @see Tools#minmax(int[],int,int)
 */
final class MinMax
{
    public final int min;
    public final int max;
    
    MinMax(int min, int max)
    {
        this.min = min;
        this.max = max;
    }

    public int range()
    {
        return max - min;
    }
}
