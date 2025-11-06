package gedi.iTiSS.utils;

public class InterquartileRange {

    /**
     * Retrieves the interquartile range of {@code array}
     * @param array The array from which the interquartile range should be calculated
     * @return The interquartile range of {@code array}
     */
    public static double interquartileRange(double[] array, boolean inplcae) {
        if (!inplcae) {
            array = array.clone();
        }
        int left = array.length / 4;
        int right = (int)(((double)array.length / (double)4) * 3.);
        return hoaresQuickselect(array, right, true) - hoaresQuickselect(array, left, true);
    }

    /**
     * Retrievs the given {@code k}'th smallest element of the given {@code array} using the Hoare's Select algorithm.
     * Also known as quickselect. (Worst O(n^2), Best O(n), Average O(n))
     * See: https://en.wikipedia.org/wiki/Quickselect
     * @param array The array to look in
     * @param k the k'th smallest element to retrieve (k=0 is the smallest)
     * @return The {@code k}'th smallest element in {@code array}
     */
    public static double hoaresQuickselect(double[] array, int k, boolean inplace) {
        if (!inplace) {
            array = array.clone();
        }
        if (k >= array.length) {
            throw new IllegalArgumentException("k needs to be smaller than the array-length. k=" + k + ", array length:" + array.length);
        }
        return select(array, 0, array.length-1, k);
    }

    private static double select(double[] array, int left, int right, int k) {
        if (left == right) {
            return array[left];
        }
        int pivotIndex = right;
        pivotIndex = partition(array, left, right, pivotIndex);
        if (k == pivotIndex) {
            return array[k];
        }
        else if (k < pivotIndex) {
            return select(array, left, pivotIndex-1, k);
        }
        else {
            return select(array, pivotIndex+1, right, k);
        }
    }

    // This is currently the Lomuto partition scheme. Should be changed to Hoare's partition scheme in the near future
    private static int partition(double[] array, int left, int right, int pivotIndex) {
        double pivotValue = array[pivotIndex];
        array[pivotIndex] = array[right];
        array[right] = pivotValue;
        int storeIndex = left;
        for (int i = left; i < right; i++) {
            if (array[i] < pivotValue) {
                double tmp = array[storeIndex];
                array[storeIndex] = array[i];
                array[i] = tmp;
                storeIndex++;
            }
        }
        double tmp = array[right];
        array[right] = array[storeIndex];
        array[storeIndex] = tmp;
        return storeIndex;
    }
}
