package gedi.util.algorithm.string.suffixarray;

import java.util.Arrays;
import java.util.Comparator;

public class NaiveSuffixArray implements ISuffixArrayBuilder {

	@Override
	public int[] buildSuffixArray(int[] input, int start, int length) {
		Integer[] re = new Integer[length];
		for (int i=0; i<re.length; i++)
			re[i] = i;
		
		class NaiveComparator implements Comparator<Integer> {
			@Override
			public int compare(Integer o1, Integer o2) {
				int len1 = length-o1;
		        int len2 = length-o2;
		        int lim = Math.min(len1, len2);

		        int k = 0;
		        while (k < lim) {
		            int c1 = input[start+o1+k];
		            int c2 = input[start+o2+k];
		            if (c1 != c2) {
		                return Integer.compare(c1,c2);
		            }
		            k++;
		        }
		        return len1 - len2;
			}
		}
		
		Arrays.sort(re,new NaiveComparator());
		System.out.println(new NaiveComparator().compare(3, 5));
		
		int[] re2 = new int[length+1];
		for (int i=0; i<re.length; i++)
			re2[i] = re[i];
		re2[re.length] = -1;
		return re2;
	}
	
	


}
