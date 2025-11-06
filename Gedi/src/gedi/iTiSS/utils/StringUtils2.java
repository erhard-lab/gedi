package gedi.iTiSS.utils;

import gedi.util.StringUtils;
import gedi.util.functions.EI;

import java.util.ArrayList;
import java.util.List;

public class StringUtils2 {
    /**
     * Returns the {@code input} string without duplicates
     *
     * @param input the string to delete the duplicates from
     * @return the {@code input} string without duplicates
     */
    public static String uniqueCharacters(String input) {
        StringBuilder unique = new StringBuilder();
        for (int i = 0; i < input.length(); i++) {
            char current = input.charAt(i);
            if (unique.toString().indexOf(current) < 0) {
                unique.append(current);
            }
        }
        return unique.toString();
    }

    /**
     * Returns all positions of {@code c} in {@code s};
     * @param s input string
     * @param c character to search for in {@code s}
     * @return an integer-array with the positions of {@code c} in {@code s}
     */
    public static int[] indicesOfChar(String s, char c) {
        List<Integer> indices = new ArrayList<>();
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == c) {
                indices.add(i);
            }
        }
        return EI.wrap(indices).toIntArray();
    }

    /**
     * Returns the first occurrence of {@code c} in {@code s}.
     * Returns -1 if {@code c} is not in {@code s}
     * @param s input
     * @param c search char
     * @return first occurrence of {@code c} in {@code s}
     */
    public static int indexOf(String s, char c) {
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == c) {
                return i;
            }
        }
        return -1;
    }

    /**
     * Creates also reversed (if to<from).
     * Both indices are inclusive!
     * @param seq
     * @param from
     * @param to
     * @return
     */
    public static char[] extract(CharSequence seq, int from, int to) {
        if (to>=from) {
            char[] re = new char[to-from+1];
            for (int i=from; i<=to; i++)
                re[i-from] = seq.charAt(i);
            return re;
        }
        char[] re = new char[from-to+1];
        for (int i=from; i>=to; i--)
            re[from-i] = seq.charAt(i);
        return re;

    }

    /**
     * Creates also reversed (if to<from).
     * Both indices are inclusive!
     * @param seq
     * @param from
     * @param to
     * @return
     */
    public static char[] extract(char[] seq, int from, int to) {
        if (to>=from) {
            char[] re = new char[to-from+1];
            for (int i=from; i<=to; i++)
                re[i-from] = seq[i];
            return re;
        }
        char[] re = new char[from-to+1];
        for (int i=from; i>=to; i--)
            re[from-i] = seq[i];
        return re;

    }

    /**
     * Creates a char array where only the characters at the given {@code indices} remain.
     * The order of the output character array is the same order as in {@code indices}.
     * Exp: {@code seq} = "abcdefg", {@code indices} = [3,2,5], result = ['d', 'c', 'f']
     * @param seq input seq
     * @param indices indices to extract from {@code seq}
     * @return character array
     */
    public static char[] extract(CharSequence seq, int[] indices) {
        char[] out = new char[indices.length];
        for (int i = 0; i < out.length; i++) {
            out[i] = seq.charAt(indices[i]);
        }
        return out;
    }

    /**
     * Creates a char array where only the characters at the given {@code indices} remain.
     * The order of the output character array is the same order as in {@code indices}.
     * Exp: {@code seq} = "abcdefg", {@code indices} = [3,2,5], result = ['d', 'c', 'f']
     * @param seq input seq
     * @param indices indices to extract from {@code seq}
     * @return character array
     */
    public static char[] extract(char[] seq, int[] indices) {
        char[] out = new char[indices.length];
        for (int i = 0; i < out.length; i++) {
            out[i] = seq[(indices[i])];
        }
        return out;
    }

    public static String concat(String glue, String... strings) {
        return StringUtils.concat(glue, strings);
    }
}
