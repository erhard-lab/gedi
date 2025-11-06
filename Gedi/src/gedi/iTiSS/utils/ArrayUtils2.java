package gedi.iTiSS.utils;

import gedi.util.ArrayUtils;

import java.util.Arrays;
import java.util.function.Function;

public class ArrayUtils2 {
    /**
     * Smooths the input array using a moving average approach with length range.
     * The output array has the same length as the input array with the first range/2 and last range/2 values
     * being 0.
     * @param input
     * @param range
     * @return
     */
    public static double[] smooth(int[] input, int range) {
        if (range > input.length) {
            throw new IllegalArgumentException("Range needs to be bigger than the array");
        }
        double[] out = new double[input.length];

        double sum = ArrayUtils.sum(input, 0, range);

        for (int i = 0; i <= input.length - range; i++) {
            int index = i + range/2;
            out[index] = sum/(double)range;
            if (i < input.length-range) {
                sum -= input[i];
                sum += input[i+range];
            }
        }

        return out;
    }

    /**
     * Smooths the input array using a moving average approach with length range.
     * The output array has the same length as the input array with the first range/2 and last range/2 values
     * being 0.
     * @param input
     * @param range
     * @return
     */
    public static double[] smooth(double[] input, int range) {
        if (range > input.length) {
            throw new IllegalArgumentException("Range needs to be bigger than the array");
        }
        double[] out = new double[input.length];

        double sum = ArrayUtils.sum(input, 0, range);

        for (int i = 0; i <= input.length - range; i++) {
            int index = i + range/2;
            out[index] = sum/(double)range;
            if (i < input.length-range) {
                sum -= input[i];
                sum += input[i+range];
            }
        }

        return out;
    }

    public static double[] subtract(double[] ary, double sub) {
        double[] out = new double[ary.length];
        for (int i = 0; i < ary.length; i++) {
            out[i] = ary[i] - sub;
        }
        return out;
    }

    public static double mean(double[] ary) {
        if (ary.length == 0) {
            return 0;
        }
        return Arrays.stream(ary).sum() / ary.length;
    }

    /**
     * Concatenates an array using {@link Object#toString()} for each item
     * and glue as separator.
     *
     * @param glue the string that is inserted between each two items
     * @param array the array
     * @return concatenation of the array
     */
    public static String concat(String glue, double[] array, int offset, int end) {
        if (array.length==0) return "";
        StringBuilder sb = new StringBuilder();
        for (int i=offset; i<end-1; i++) {
            sb.append(array[i]);
            sb.append(glue);
        }
        sb.append(array[end-1]);
        return sb.toString();
    }

    /**
     * Concatenates an array using {@link Object#toString()} for each item
     * and glue as separator.
     *
     * @param glue the string that is inserted between each two items
     * @param array the array
     * @return concatenation of the array
     */
    public static String concat(String glue, double[] array) {
        return concat(glue, array, 0, array.length);
    }

    public static double[] subtract(double[] total, double[] sub) {
        if (total==null||total.length!=sub.length) total = new double[sub.length];
        for (int i=0; i<total.length; i++)
            total[i]-=sub[i];
        return total;
    }

    public static double[] divide(double[] a, double b) {
        double[] re = new double[a.length];
        for (int i=0; i<a.length; i++)
            re[i] = a[i]/b;
        return re;
    }

    public static boolean contains(int[] ary, int x) {
        for (int i : ary) {
            if (i == x) {
                return true;
            }
        }
        return false;
    }

    public static boolean contains(String[] ary, String x) {
        for (String s : ary) {
            if (s.equals(x)) {
                return true;
            }
        }
        return false;
    }
}
