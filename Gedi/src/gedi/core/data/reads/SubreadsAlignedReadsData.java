package gedi.core.data.reads;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import gedi.app.Gedi;
import gedi.core.data.reads.AlignedReadsDataFactory.VarIndel;
import gedi.core.data.reads.subreads.PairedEndToSubreadsConverter;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.sparse.AutoSparseDenseIntArrayCollector;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.BinaryReader;

/**
 * This is a data object once the raw data was summarized (i.e. if is more than one read for the same RNA fragment because of paired-end or UMIs)
 * Summarized means that mismatches have already been made consistent (and moved to the sense strand), and there is no geometry information anymore.
 * For these, the mapped read is composed of subreads (e.g. for paired end sequencing, there are two or three subreads, sense and antisense, or sense, overlap and antisense)
 * This also means that if sequencing was strand unspecific, the read mapping has been moved to the sense strand!
 * 
 * The exact semantics of the subreads (as the above mentioned sense, antisense and overlap) is  defined in a separate object implementing {@link SubreadsSemantic}
 * 
 * This influences on distinct sequences: the subread ids are distinct sequence specific,i.e. there can be two distinct subread patterns matching to the same genomic region; 
 * if two reads with the same mismatch pattern have a distinct subread pattern, we have two distinct sequences (this can happen, geometries are distinct sequences specific).
 * 
 * For paired end, this means: Two distinct sequence might be merged if and only if they have the same geometry, and the mismatches that make them distinct are pair inconsistent and therefore removed. This is done via the fac.makeDistinct() call int {@link PairedEndToSubreadsConverter}
 * however, it is impossible that a new distinct sequence is created. 
 * 
 * For 10x, this is different: In most of the cases, the reads for one UMI will be very different to the reads of a different UMI (even if the UMI covers the same region and has the same mismatch pattern). Thus, it is to be expected that each distinct sequence will only have a single count!
 * 
 * @author flo
 *
 */
public class SubreadsAlignedReadsData extends DefaultAlignedReadsData implements HasSubreads {

	// semindex, end of subread, semindex, end of subread, semindex, ..., end of subread, semindex
	int[][] subreadGeom;
	int[][] gaps;

	public SubreadsAlignedReadsData() {}


	public int getNumSubreads(int distinct) {
		return subreadGeom[distinct].length/2+1;
	}
	
	public int getSubreadStart(int distinct, int index) {
		if (index==0) return 0;
		return subreadGeom[distinct][index*2-1];
	}
	public int getSubreadId(int distinct, int index) {
		return subreadGeom[distinct][index*2];
	}
	@Override
	public IntIterator getGapPositions(int distinct) {
		return EI.wrap(gaps[distinct]);
	}

	
	@Override
	public HasSubreads asSubreads(boolean sense) {
		return this;
	}
	
	
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		super.deserialize(in);
		int d = getDistinctSequences();
		
		subreadGeom = new int[d][];
		gaps = new int[d][];
		for (int i=0; i<d; i++) {
			int numsub = in.getCInt();
			subreadGeom[i] = new int[numsub*2-1];
			for (int s=0; s<subreadGeom[i].length; s++) 
					subreadGeom[i][s]=in.getCInt();
			int numgap = in.getCInt();
			gaps[i] = new int[numgap];
			for (int s=0; s<gaps[i].length; s++) 
				gaps[i][s]=in.getCInt();
		}
		
	}
	
//	@Override
//	public String toString() {
//		StringBuilder sb = new StringBuilder();
//		for (int d=0; d<getDistinctSequences(); d++) {
//			if (hasId()) 
//				sb.append(getId(d)).append(": ");
//			if (hasNonzeroInformation()) {
//				sb.append("[");
//				int[] inds = getNonzeroCountIndicesForDistinct(d);
//				for (int i=0; i<inds.length; i++)
//					sb.append(inds[i]+":"+getNonzeroCountValueForDistinct(d, i)+",");
//				sb.deleteCharAt(sb.length()-1);
//				sb.append("]");
//			} else
//				sb.append(Arrays.toString(getCountsForDistinctInt(d, ReadCountMode.All)));
//			sb.append(" x");
//			sb.append(getMultiplicity(d));
//			if (hasWeights()) 
//				sb.append(" (w=").append(String.format("%.2f", getWeight(d))).append(")");
//			sb.append(" ");
//			addSubreadToString(sb, d);
//			for (AlignedReadsVariation var : getVariations(d))
//				sb.append("\t"+var);
//			
//			if (d<getDistinctSequences()-1) sb.append(" ~ ");
//		}
//		return sb.toString();
//	}


	@Override
	public boolean isConsistentlyContained(ReferenceGenomicRegion<?> read, ReferenceGenomicRegion<?> reference, int d) {
		return HasSubreads.super.isConsistentlyContained(read, reference, d);
	}
	@Override
	public boolean isFalseIntron(int pos, int distinct) {
		return HasSubreads.super.isFalseIntron(pos, distinct);
	}
	@Override
	public int getNumParts(ReferenceGenomicRegion<?> read, int distinct) {
		return HasSubreads.super.getNumParts(read, distinct);
	}
	
	/**
	 * vars must be sorted
	 * @param numCond
	 * @param condition
	 * @param subread
	 * @param gap
	 * @param vars
	 * @return
	 */
	public static SubreadsAlignedReadsData createUnitCount(int numCond, int condition, int[] subread, int[] gap, List<VarIndel> vars) {
		SubreadsAlignedReadsData re = new SubreadsAlignedReadsData();
		
		re.conditions = numCond;
		re.count = new int[1][1]; re.count[0][0]=1;
		re.nonzeros = new int[1][1]; re.nonzeros[0][0]=condition;
		re.var = new int[1][vars.size()];
		for (int i=0; i<vars.size(); i++) 
			re.var[0][i] = vars.get(i).var;
		re.indels = new CharSequence[1][vars.size()];
		for (int i=0; i<vars.size(); i++)
			re.indels[0][i] = vars.get(i).indel;
		re.multiplicity = new int[] {1};
		re.subreadGeom = new int[][] {subread};
		re.gaps = new int[][] {gap};

		return re;
	}
	
	public static SubreadsAlignedReadsData create(int numCond, int[][] count, int[][] nonzeros, int[][] vars, CharSequence[][] indels, int[][] subreadGeom, int[][] gaps) {
		SubreadsAlignedReadsData re = new SubreadsAlignedReadsData();
		
		re.conditions = numCond;
		re.count = count;
		re.nonzeros = nonzeros;
		re.var = vars;
		re.indels = indels;
		re.multiplicity = new int[count.length];
		Arrays.fill(re.multiplicity, 1);
		re.subreadGeom = subreadGeom;
		re.gaps = gaps;

		return re;
	}
	
	/**
	 * 
	 * @param mapping is oldconditions->newconditions
	 * @return
	 */
	public SubreadsAlignedReadsData selectMergeConditions(int numNewConditions, int[][] mapping) {
		int[][] nonzerore = null;
		int[][] countre = new int[getDistinctSequences()][];
		
		if (hasNonzeroInformation()) {
			if (numNewConditions>5) {
				nonzerore = new int[getDistinctSequences()][];
				IntArrayList nonzerol = new IntArrayList(3);
				IntArrayList countl = new IntArrayList(3);
				for (int d=0; d<getDistinctSequences(); d++)  {
					for (int n=0; n<nonzeros[d].length; n++) {
						int oldcond = nonzeros[d][n];
						int count = this.count[d][n];
						
						for (int newcond : mapping[oldcond]) {
							nonzerol.add(newcond);
							countl.add(count);
						}
					}
					
					// now we collected everything, but it might not be unique (and is not sorted)
					int[] nonzeroraw = nonzerol.getRaw();
					int[] countraw = countl.getRaw();
					ArrayUtils.parallelSort(nonzeroraw, countraw,0,nonzerol.size());
					int newlen = makeUnique(nonzerol.size(), nonzeroraw,countraw);
					
					nonzerore[d] = ArrayUtils.slice(nonzeroraw, 0, newlen);
					countre[d] = ArrayUtils.slice(countraw, 0, newlen);
					
					countl.clear();
					nonzerol.clear();
				}
			} else {
				for (int d=0; d<getDistinctSequences(); d++)  {
					countre[d] = new int[numNewConditions];
					
					for (int n=0; n<nonzeros[d].length; n++) {
						int oldcond = nonzeros[d][n];
						int count = this.count[d][n];
						
						for (int newcond : mapping[oldcond]) {
							countre[d][newcond]+=count;							
						}
					}
				}
			}
		} else {
			for (int d=0; d<getDistinctSequences(); d++)  {
				countre[d] = new int[numNewConditions];
				
				for (int n=0; n<count.length; n++) {
					int oldcond = n;
					int count = this.count[d][n];
					
					for (int newcond : mapping[oldcond]) {
						countre[d][newcond]+=count;							
					}
				}
			}
		}
		
		return SubreadsAlignedReadsData.create(numNewConditions, countre, nonzerore, this.var, this.indels, this.subreadGeom, this.gaps);
	}
	
	private static int makeUnique(int len, int[] s, int[] c) {
		int index = 0;
		for (int i=index+1; i<len; i++) {
			if (s[i]!=s[index]) {
				index++;
				s[index]=s[i];
				c[index]=c[i];
			} else {
				c[index]+=c[i];
			}
		}
		return index+1;
	}


	public static int compare(SubreadsAlignedReadsData a, int ai, SubreadsAlignedReadsData b, int bi) {
		int re = ArrayUtils.compare(a.var[ai],b.var[bi]);
		if (re!=0) return re;
		return ArrayUtils.compare(a.subreadGeom[ai], b.subreadGeom[bi]);
	}


	public static void main(String[] args) throws IOException {
		Gedi.startup(false);
		GenomicRegionStorage<DefaultAlignedReadsData> st = Gedi.load(args[0]);
		PairedEndToSubreadsConverter conv = new PairedEndToSubreadsConverter(true);
		
		for (ImmutableReferenceGenomicRegion<DefaultAlignedReadsData> r : st.ei().loop()) {
			System.out.println(r);
			System.out.println(r.getData().getNumParts(r, 0));
			System.out.println(conv.convert(r, true,null,null));
			System.out.println(conv.convert(r, true,null,null).getData().getNumParts(conv.convert(r, true,null,null), 0));
			System.out.println(conv.convert(r, false,null,null));
			System.out.println(conv.convert(r, false,null,null).getData().getNumParts(conv.convert(r, false,null,null), 0));
			System.out.println();
			
		}
		
		
	}


	public static void setVarsAndSubreads(SubreadsAlignedReadsData d, int inIndex, int[][] vars,
			CharSequence[][] indels, int[][] subreadGeom, int[][] gaps, int outIndex) {
		vars[outIndex] = d.var[inIndex];
		indels[outIndex] = d.indels[inIndex];
		subreadGeom[outIndex] = d.subreadGeom[inIndex];
		gaps[outIndex] = d.gaps[inIndex];
	}

	
	

}
