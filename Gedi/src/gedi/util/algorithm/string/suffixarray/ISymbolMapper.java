package gedi.util.algorithm.string.suffixarray;

/**
 * Symbol mappers (reversible int-coding).
 */
interface ISymbolMapper
{
    void map(int [] input, int start, int length);
    void undo(int [] input, int start, int length);
}
