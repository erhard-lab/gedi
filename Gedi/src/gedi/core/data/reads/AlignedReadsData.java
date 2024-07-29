package gedi.core.data.reads;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;

import gedi.app.extension.GlobalInfoProvider;
import gedi.core.data.HasConditions;
import gedi.core.data.reads.AlignedReadsDataFactory.VarIndel;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.FileUtils;
import gedi.util.FunctorUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.datastructure.array.sparse.AutoSparseDenseIntArrayCollector;
import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.IntDoubleToDoubleFunction;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;
import gedi.util.math.stat.Ranking;
import gedi.util.sequence.DnaSequence;


/**
 * <p>
 * Usually used as data for {@link GenomicRegion}. The GenomicRegion should indicate where a read originated, i.e. gaps due to biological reasons
 * like introns should be modelled using an intron in the GenomicRegion. The Variations here should be used to indicate gaps and mismatches due
 * to technical reasons.
 * </p>
 * 
 * <p>
 * Positions refer to the covered part of the genome sequence (not the alignment nor the read sequence) in 5' to 3' direction.
 * 
 * <p><pre>
 * Example:
 * Genome            AAAAAAAAAATTTTTTTTAAATTTTTTTTTTT-CCCCCCCCCCCCCAAAAAAAAAAAATTTTTTTTTTTTTTACCCCCCCCCCCCCCCCAAAAAAAAAAAAAA
 * Read 1 alignment            TTTTTTTT---TTTTTTTTTTTACCCCCCCCCCCCC------------TTTTTTTTTTTTTTGCCCCCCCCCCCCCCCC
 * Read 2 alignment            TTTTTTTT---TTTTTTTTTTT-CCCCCCCCCCCCC------------TTTTTTTTTTTTTTACCCCCCCCCCCCCCCC
 * Read 3 alignment            TTTTTTTTAAATTTTTTTTTTT-CCCCCCCCCCCCC------------TTTTTTTTTTTTTTACCCCCCCCCCCCCCCC
 * </pre></p>
 * 
 * <p>
 * The long gap in the read is an intron, so the {@link GenomicRegion} of the read should be 10-46|58-81 (assuming plus strand).
 * Here, we have 3 distinct sequences ({@link #getDistinctSequences()} is 3). 
 * The first has 3 variations ({@link #getVariationCount(int)} is 3 for index 0). The first is a deletion (isDeletion(0,0) is true) at position
 * 8 (getDeletionPos(0,0) is 8) with length 3 (getDeletionLength(0,0) is 3). The second is an insertion
 * at position 22 with sequence A  and the third is a mismatch at position 49 with genomic char A and read char G (positions refer to covered parts of the genome)
 * The second read has one variation (the deletion), the last none.
 * </p>
 * 
 * <p>
 * If reads are aligned to the negative genomic strand, everything is considered in 5' to 3' direction. I.e. in genomic direction, variation indices
 * are right-to-left and positions are alignmentendposition-genomic position (introns must also be considered)
 * </p>
 * 
 * <p>12/9/15: Each distinct sequence may have an integer id. If read from BinaryFile, having or not having an id is decided based on an attribte in the reader's context.
 * The returned id is -1 otherwise.
 * </p>
 * <p>31/1/18: second read: the genomic base is the one from the strand of the second read (i.e. from the opposite strand of the rgr); all read sequences naturally refer to the sequence of the read itself (also opposite strand of rgr, so to say); softclips are also defined by the direction, i.e. 3p softclip is a the end of the second reads (i.e. in 5' direction of rgr)
 * </p>
 * 
 * Geometry information is also w.r.t. 5' to 3'. I.e. before on - is the high position
 * 
 * If N is sequenced, a A->A mismatch is stored (if the genomic base is an A)!
 * If N is sequenced and softclipped, an A is stored!
 * 
 * @author erhard
 *
 */
public interface AlignedReadsData extends BinarySerializable, GlobalInfoProvider, HasConditions {

	public static final String HASIDATTRIBUTE = "HASIDATTRIBUTE";
	public static final String HASWEIGHTATTRIBUTE = "HASWEIGHT";
	public static final String HASGEOMETRYATTRIBUTE = "HASGEOMETRY";
	public static final String CONDITIONSATTRIBUTE = "CONDITIONS";
	public static final String SPARSEATTRIBUTE = "SPARSE";
	
	public static final String BETTERPOSATTRIBUTE = "BPOS";
	
	
	public static final String BARCODEATTRIBUTE = "BARCODE";
	
	
	default DynamicObject getGlobalInfo() {
		LinkedHashMap<String, Integer> map = new LinkedHashMap<String, Integer>();
		map.put(HASIDATTRIBUTE, hasId()?1:0);
		map.put(HASWEIGHTATTRIBUTE, hasWeights()?1:0);
		map.put(HASGEOMETRYATTRIBUTE, hasGeometry()?1:0);
		map.put(CONDITIONSATTRIBUTE, getNumConditions());
		map.put(SPARSEATTRIBUTE, getNumConditions()>5?2:0); // barcoded reads cannot be handled in a sparse manner by the factory, but they can be saved sparsely!
//		map.put(SPARSEATTRIBUTE, hasNonzeroInformation()?2:0);
		if (this instanceof HasBarcodes)
			map.put(BARCODEATTRIBUTE, ((HasBarcodes)this).getBarcodeLength());
		
		map.put(BETTERPOSATTRIBUTE, 1);
		return DynamicObject.from(map);
	}
	
	/**
	 * If this read is only observed w/o mismatches, 1 is returned.
	 * @return
	 */
	int getDistinctSequences();
	int getCount(int distinct, int condition);
	int getVariationCount(int distinct); 
	
	int getId(int distinct);

	float getWeight(int distinct);
	boolean hasWeights();
	boolean hasGeometry();
	int getGeometryBeforeOverlap(int distinct);
	int getGeometryOverlap(int distinct);
	int getGeometryAfterOverlap(int distinct);
	int getRawGeometry(int distinct);
	
	
	default int[] getNonzeroCountIndicesForDistinct(int distinct) {
		return null;
	}
	
	default int getNonzeroCountValueForDistinct(int distinct, int index) {
		return -1;
	}
	
	default double getNonzeroCountValue(int distinct, int index, ReadCountMode mode) {
		return mode.computeCount(getNonzeroCountValueForDistinct(distinct,index), getMultiplicity(distinct), getWeight(distinct));
	}
	
	default int getNonzeroCountValueInt(int distinct, int index, ReadCountMode mode) {
		return mode.computeCountInt(getNonzeroCountValueForDistinct(distinct,index), getMultiplicity(distinct), getWeight(distinct));
	}
	default int getNonzeroCountValueFloor(int distinct, int index, ReadCountMode mode) {
		return mode.computeCountFloor(getNonzeroCountValueForDistinct(distinct,index), getMultiplicity(distinct), getWeight(distinct));
	}
	
	default boolean hasNonzeroInformation() {
		return false;
	}
	
	default int getMappedLength(int distinct) {
		return getGeometryBeforeOverlap(distinct)+getGeometryOverlap(distinct)+getGeometryAfterOverlap(distinct);
	}
	
	default int getMappedLengthRead1(int distinct) {
		return getGeometryBeforeOverlap(distinct)+getGeometryOverlap(distinct);
	}
	
	default int getMappedLengthRead2(int distinct) {
		return getGeometryAfterOverlap(distinct)+getGeometryOverlap(distinct);
	}
	
	boolean isVariationFromSecondRead(int distinct, int index);
	
	default int getFirstReadClip(int d) {
		int re = 0;
		for (int v=0; v<getVariationCount(d); v++)
			if (!isVariationFromSecondRead(d, v) && isSoftclip(d, v))
				re+=getSoftclip(d, v).length();
		return re;
	}
	
	default int getSecondReadClip(int d) {
		int re = 0;
		for (int v=0; v<getVariationCount(d); v++)
			if (isVariationFromSecondRead(d, v) && isSoftclip(d, v))
				re+=getSoftclip(d, v).length();
		return re;
	}
	
	
	default boolean hasMismatch(int distinct) {
		for (int v=0; v<getVariationCount(distinct); v++)
			if (isMismatch(distinct, v))
				return true;
		return false;
	}
	
	default boolean hasMismatch(int distinct, char genomic, char read) {
		for (int v=0; v<getVariationCount(distinct); v++)
			if (isMismatch(distinct, v) && getMismatchGenomic(distinct, v).charAt(0)==genomic && getMismatchRead(distinct, v).charAt(0)==read)
				return true;
		return false;
	}
	
	/**
	 * Does not check strand!
	 * @param read
	 * @param reference
	 * @param d
	 * @return
	 */
	default boolean isConsistentlyContained(ReferenceGenomicRegion<?> read,
			ReferenceGenomicRegion<?> reference, int d) {
		if (!hasGeometry())
			return reference.getRegion().containsUnspliced(read.getRegion());
		
		return reference.getRegion().containsUnspliced(extractRead1(read, d).getRegion())
				&& reference.getRegion().containsUnspliced(extractRead2(read, d).getRegion());
	}

	
	default ReferenceGenomicRegion<Void> extractRead1(ReferenceGenomicRegion<?> read,int d) {
		return new ImmutableReferenceGenomicRegion<>(
				read.getReference(), 
				read.map(new ArrayGenomicRegion(0,getMappedLengthRead1(d)))
				);
	}
	
	
	default ReferenceGenomicRegion<Void> extractRead2(ReferenceGenomicRegion<?> read,int d) {
		return new ImmutableReferenceGenomicRegion<>(
				read.getReference(), 
				read.map(new ArrayGenomicRegion(read.getRegion().getTotalLength()-getMappedLengthRead2(d),read.getRegion().getTotalLength()))
				);
	}
	
	default int getNumParts(ReferenceGenomicRegion<?> read,int d) {
		if (!hasGeometry() || getGeometryOverlap(d)>0) return read.getRegion().getNumParts();
		int p1 = read.map(getGeometryBeforeOverlap(d)-1);
		int p2 = read.map(getGeometryBeforeOverlap(d));
		if (Math.abs(p1-p2)!=1)
			return read.getRegion().getNumParts()-1;
		return read.getRegion().getNumParts();
	}
	
	/**
	 * l must be in mapped region coordinates
	 * @param l
	 * @param d
	 * @return
	 */
	default boolean isFalseIntron(int l, int d) {
		return hasGeometry() && l==getGeometryBeforeOverlap(d);
	}

	
	default int getReadLength(int d, int regionLength) {
		int l = regionLength;
		for (int v=0; v<getVariationCount(d); v++)
			if (isInsertion(d, v))
				l+=getInsertion(d, v).length();
			else if (isDeletion(d, v))
				l-=getDeletion(d, v).length();
			else if (isSoftclip(d, v))
				l+=getSoftclip(d, v).length();
		
		return l;
	}
	
	/**
	 * Gets the end position of the first read in the mapped region
	 * @param d
	 * @param pos
	 * @param readLength1
	 * @return
	 */
	default int getReadLength1(int d) {
		int l = getMappedLengthRead1(d);
		for (int v=0; v<getVariationCount(d); v++)
			if (isInsertion(d, v) && !isVariationFromSecondRead(d, v))
				l+=getInsertion(d, v).length();
			else if (isDeletion(d, v) && !isVariationFromSecondRead(d, v))
				l-=getDeletion(d, v).length();
			else if (isSoftclip(d, v) && !isVariationFromSecondRead(d, v))
				l+=getSoftclip(d, v).length();
		
		return l;
	}
	
	/**
	 * Gets the length of the second read in the mapped region
	 * @param d
	 * @param pos
	 * @param readLength2
	 * @return
	 */
	default int getReadLength2(int d) {
		int l = getMappedLengthRead2(d);
		for (int v=0; v<getVariationCount(d); v++)
			if (isInsertion(d, v) && isVariationFromSecondRead(d, v))
				l+=getInsertion(d, v).length();
			else if (isDeletion(d, v) && isVariationFromSecondRead(d, v))
				l-=getDeletion(d, v).length();
			else if (isSoftclip(d, v) && isVariationFromSecondRead(d, v))
				l+=getSoftclip(d, v).length();
		
		return l;
	}
	
	

	/**
	 * Maps from the coordinate system of the read mapping to read coordinates. Read coordinates is as if both read sequences had been concatenated (second reversed!)!
	 * (paying attention to any insertion, softclips and mapping to genomic coordinates)
	 * <br>
	 * Caution: Is pos is an {@link #getInsertionPos(int, int)}, the end of the insertion is reported! This is necessary, because
	 * an insertion and mismatch can have the same position (the mismatch occurrs right after the insertion)
	 * 
	 * @param d
	 * @param pos
	 * @param secondReadInOverlap if position is in overlap, report second read position?
	 * @return
	 */
	default int mapToRead(int d, int pos, boolean secondReadInOverlap) {
		return mapToRead(d, pos, secondReadInOverlap,-1);
	}
	
	/**
	 * Maps from the coordinate system of the read mapping to read coordinates. Read coordinates is as if both read sequences had been concatenated (second reversed!)!
	 * (paying attention to any insertion, softclips and mapping to genomic coordinates)
	 * <br>
	 * Caution: Is pos is an {@link #getInsertionPos(int, int)}, the end of the insertion is reported! This is necessary, because
	 * an insertion and mismatch can have the same position (the mismatch occurrs right after the insertion)
	 * 
	 * <br>
	 * Here the read lengths should be given (important, if reads have been hardclipped (e.g. adapter trimmed) to have the first position of the second read always at the same position)
	 * 
	 * @param d
	 * @param pos
	 * @param secondReadInOverlap if position is in overlap, report second read position?
	 * @return
	 */
	default int mapToRead(int d, int pos, boolean secondReadInOverlap, int readLength12) {
		boolean inFirst = isPositionInFirstRead(d, pos);
		boolean inSecond = isPositionInSecondRead(d, pos);
		
		boolean reportFirstRead = inFirst && (!inSecond || !secondReadInOverlap);
		
		if (reportFirstRead) 
			return mapToRead1(d, pos);
		
		return mapToRead2(d, pos, readLength12);
	}
	
	default boolean isPositionInFirstRead(int d, int pos) {
		if (!hasGeometry()) return true;
		int first = getMappedLengthRead1(d);
		return pos<first;
	}
	
	default boolean isPositionInOverlap(int d, int pos) {
		if (!hasGeometry()) return false;
		return isPositionInFirstRead(d, pos) && isPositionInSecondRead(d, pos);
	}
	
	default boolean isPositionInSecondRead(int d, int pos) {
		if (!hasGeometry()) return false;
		int second = getMappedLengthRead2(d);
		return pos>=getMappedLength(d)-second;
	}
	
	/**
	 * Returns -1 if the position is in a deletion!
	 * @param d
	 * @param pos
	 * @return
	 */
	default int mapToRead1(int d, int pos) {
		int add = 0;
		for (int v=0; v<getVariationCount(d); v++)
			if (isInsertion(d, v) && !isVariationFromSecondRead(d, v) && getInsertionPos(d, v)<=pos)
				add+=getInsertion(d, v).length();
			else if (isDeletion(d, v) && !isVariationFromSecondRead(d, v) && getDeletionPos(d, v)<=pos) {
				if (pos<getDeletionPos(d, v)+getDeletion(d, v).length())
					return -1;
				add-=getDeletion(d, v).length();
			} else if (isSoftclip(d, v) && isSoftclip5p(d, v) && !isVariationFromSecondRead(d, v))
				add+=getSoftclip(d, v).length();
		return pos+add;
	}
	/**
	 * Returns -1 if the position is in a deletion!
	 * @param d
	 * @param pos
	 * @return
	 */
	default int mapToRead2(int d, int pos) {
		return mapToRead2(d, pos, -1);
	}
	/**
	 * Returns -1 if the position is in a deletion!
	 * @param d
	 * @param pos
	 * @return
	 */
	default int mapToRead2(int d, int pos, int readLength12) {
		if (readLength12<0) readLength12 = getReadLength1(d)+getReadLength2(d);
		
		int add = readLength12-getMappedLength(d);
		
		for (int v=getVariationCount(d)-1; v>=0; v--)
			if (isInsertion(d, v) && isVariationFromSecondRead(d, v) && getInsertionPos(d, v)>pos)
				add-=getInsertion(d, v).length();
			else if (isDeletion(d, v) && isVariationFromSecondRead(d, v) && getDeletionPos(d, v)+getDeletion(d, v).length()>pos) {
				if (pos>=getDeletionPos(d, v))
					return -1;
				add+=getDeletion(d, v).length();
			}
			else if (isSoftclip(d, v) && !isSoftclip5p(d, v) && isVariationFromSecondRead(d, v))
				add-=getSoftclip(d, v).length();
		return pos+add;
	}
	



	default boolean isMismatchN(int distinct, int index) {
		return StringUtils.equals(getMismatchGenomic(distinct, index),getMismatchRead(distinct, index));
	}
	boolean isMismatch(int distinct, int index);
	int getMismatchPos(int distinct, int index);
	CharSequence getMismatchGenomic(int distinct, int index);
	CharSequence getMismatchRead(int distinct, int index);

	boolean isSoftclip(int distinct, int index);
	boolean isSoftclip5p(int distinct, int index);
	/**
	 * This is a read sequence
	 * @param distinct
	 * @param index
	 * @return
	 */
	CharSequence getSoftclip(int distinct, int index);

	
	/**
	 * Gap in reference sequence
	 * @param distinct
	 * @param index
	 * @return
	 */
	boolean isInsertion(int distinct, int index);
	int getInsertionPos(int distinct, int index);
	
	/**
	 * This is a read sequence
	 * @param distinct
	 * @param index
	 * @return
	 */
	CharSequence getInsertion(int distinct, int index);

	boolean isDeletion(int distinct, int index);
	int getDeletionPos(int distinct, int index);
	/**
	 * This is a genomic sequence!
	 * @param distinct
	 * @param index
	 * @return
	 */
	CharSequence getDeletion(int distinct, int index);


	/**
	 * 0 means unknown
	 * @param distinct
	 * @return
	 */
	int getMultiplicity(int distinct);

	
	
	default boolean hasId() {
		return getId(0)>=0;
	}
	
	default boolean isAmbigousMapping(int distinct) {
		return getMultiplicity(distinct)!=1;
	}
	
	default boolean isUniqueMapping(int distinct) {
		return !isAmbigousMapping(distinct);
	}
	
	default int positionToGenomic(int pos, ReferenceSequence ref, GenomicRegion region) {
		if (ref.getStrand()==Strand.Minus)
			return region.map(region.getTotalLength()-1-pos);
		return region.map(pos);
	}
	
	
	default double getCount(int distinct, int condition, ReadCountMode mode) {
		return mode.computeCount(getCount(distinct, condition), getMultiplicity(distinct), getWeight(distinct));
	}
	
	
	
	default double getTotalCountForCondition(int condition, ReadCountMode mode) {
		double re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getCount(i, condition, mode);
		return re;
	}
	
	default double[] getCountsForDistinct(int distinct, ReadCountMode mode) {
		double[] re = new double[getNumConditions()];
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] = getNonzeroCountValue(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] = getCount(distinct, i, mode);
		
		return re;
	}
	
	default int[] getCountsForDistinctInt(int distinct, ReadCountMode mode) {
		int[] re = new int[getNumConditions()];
		
		if (hasNonzeroInformation())  {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] = getNonzeroCountValueInt(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] = getCountInt(distinct, i, mode);
		return re;
	}
	
	default int[] getCountsForDistinctFloor(int distinct, ReadCountMode mode) {
		int[] re = new int[getNumConditions()];
		
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] = getNonzeroCountValueFloor(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] = getCountFloor(distinct, i, mode);
		return re;
	}
	
	
	default double[] addCountsForDistinct(int distinct, double[] re, ReadCountMode mode) {
		if (re==null)
			re = new double[getNumConditions()];
		else if (re.length!=getNumConditions())
			throw new RuntimeException("Length does not match!");
	
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] += getNonzeroCountValue(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] += getCount(distinct, i, mode);
		return re;
	}

	/**
	 * do not add as is to re, but at condition index i of this at position reindex[i] in re! 
	 * @param distinct
	 * @param re
	 * @param reindex
	 * @param mode
	 * @return
	 */
	default double[] addCountsForDistinct(int distinct, double[] re, int[] reindex, ReadCountMode mode) {
		if (re==null)
			throw new RuntimeException("Null not allowed!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[reindex[inds[i]]] += getNonzeroCountValue(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[reindex[i]] += getCount(distinct, i, mode);
		return re;
	}

	
	default NumericArray addCountsForDistinct(int distinct, NumericArray re, ReadCountMode mode) {
		
		if (re==null)
			re = NumericArray.createMemory(getNumConditions(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getNumConditions() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				mode.addCount(re,inds[i],getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int c=0; c<re.length(); c++)
				mode.addCount(re,c,getCount(distinct, c),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}
	
	default AutoSparseDenseDoubleArrayCollector addCountsForDistinct(int distinct, AutoSparseDenseDoubleArrayCollector re, ReadCountMode mode) {
		
		if (re==null)
			re = new AutoSparseDenseDoubleArrayCollector(Math.max(50, getNumConditions()>>3),getNumConditions());
		else if (re.length()!=getNumConditions())
			throw new RuntimeException("Length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				mode.addCount(re,inds[i],getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int c=0; c<re.length(); c++)
				mode.addCount(re,c,getCount(distinct, c),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}

	/**
	 * do not add as is to re, but for condition index i of this at position reindex[i] in re! 
	 * @param distinct
	 * @param re
	 * @param reindex
	 * @param mode
	 * @return
	 */
	default AutoSparseDenseDoubleArrayCollector addCountsForDistinct(int distinct, AutoSparseDenseDoubleArrayCollector re, int[] reindex, ReadCountMode mode) {
		
		if (re==null)
			throw new RuntimeException("Null not allowed!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				mode.addCount(re,reindex[inds[i]],getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int c=0; c<re.length(); c++)
				mode.addCount(re,reindex[c],getCount(distinct, c),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}

	/**
	 * do not add as is to re, but, for condition index i of this at each position reindex[i][.] in re! 
	 * @param distinct
	 * @param re
	 * @param reindex
	 * @param mode
	 * @return
	 */
	default AutoSparseDenseDoubleArrayCollector addCountsForDistinct(int distinct, AutoSparseDenseDoubleArrayCollector re, int[][] reindex, ReadCountMode mode) {
		
		if (re==null)
			throw new RuntimeException("Null not allowed!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				for (int toind : reindex[inds[i]])
					mode.addCount(re,toind,getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int c=0; c<re.length(); c++)
				for (int toind : reindex[c])
					mode.addCount(re,toind,getCount(distinct, c),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}

	default AutoSparseDenseIntArrayCollector addCountsForDistinctInt(int distinct, AutoSparseDenseIntArrayCollector re, ReadCountMode mode) {
		
		if (re==null)
			re = new AutoSparseDenseIntArrayCollector(Math.max(50, getNumConditions()>>3),getNumConditions());
		else if (re.length()!=getNumConditions())
			throw new RuntimeException("Length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				mode.addCountInt(re,inds[i],getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int c=0; c<re.length(); c++)
				mode.addCountInt(re,c,getCount(distinct, c),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}

	
	default int[] addCountsForDistinctInt(int distinct, int[] re, ReadCountMode mode) {
		if (re==null)
			re = new int[getNumConditions()];
		else if (re.length!=getNumConditions())
			throw new RuntimeException("Length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] += getNonzeroCountValueInt(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] += getCountInt(distinct, i, mode);
		return re;
	}
	
	/**
	 * do not add as is to re, but at condition index i of this at position reindex[i] in re! 
	 * @param distinct
	 * @param re
	 * @param reindex
	 * @param mode
	 * @return
	 */
	default long[] addCountsForDistinctInt(int distinct, long[] re, int[] reindex, ReadCountMode mode) {
		if (re==null)
			throw new RuntimeException("Null not allowed!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[reindex[inds[i]]] += getNonzeroCountValue(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++) {
				re[reindex[i]] += getCount(distinct, i, mode);
				System.out.println(i+"->"+reindex[i]);
			}
		return re;
	}

	
	default long[] addCountsForDistinctInt(int distinct, long[] re, ReadCountMode mode) {
		if (re==null)
			re = new long[getNumConditions()];
		else if (re.length!=getNumConditions())
			throw new RuntimeException("Length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] += getNonzeroCountValueInt(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] += getCountInt(distinct, i, mode);
		return re;
	}
	
	default int[] addCountsForDistinctFloor(int distinct, int[] re, ReadCountMode mode) {
		if (re==null)
			re = new int[getNumConditions()];
		else if (re.length!=getNumConditions())
			throw new RuntimeException("Length does not match!");
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re[inds[i]] += getNonzeroCountValueFloor(distinct, i, mode);
		}
		else
			for (int i=0; i<re.length; i++)
				re[i] += getCountFloor(distinct, i, mode);
		return re;
	}
	
	default NumericArray getCountsForDistinct(NumericArray re, int distinct, ReadCountMode mode) {
		if (re==null)
			re = NumericArray.createMemory(getNumConditions(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getNumConditions() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");
		
		
		if (hasNonzeroInformation()) {
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				mode.getCount(re,inds[i],getNonzeroCountValueForDistinct(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		}
		else
			for (int i=0; i<re.length(); i++)
				mode.getCount(re,i,getCount(distinct, i),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}
	
	
	default double[] addCountsForCondition(int condition, double[] re, ReadCountMode mode) {
		if (re==null)
			re = new double[getDistinctSequences()];
		else if (re.length!=getDistinctSequences())
			throw new RuntimeException("Length does not match!");
		
		for (int i=0; i<re.length; i++)
			re[i] += getCount(i, condition, mode);
		return re;
	}
	
	default int[] addCountsForConditionInt(int condition, int[] re, ReadCountMode mode) {
		if (re==null)
			re = new int[getDistinctSequences()];
		else if (re.length!=getDistinctSequences())
			throw new RuntimeException("Length does not match!");

		for (int i=0; i<re.length; i++)
			re[i] += getCountInt(i, condition, mode);
		return re;
	}
	
	default int[] addCountsForConditionFloor(int condition, int[] re, ReadCountMode mode) {
		if (re==null)
			re = new int[getDistinctSequences()];
		else if (re.length!=getDistinctSequences())
			throw new RuntimeException("Length does not match!");

		for (int i=0; i<re.length; i++)
			re[i] += getCountFloor(i, condition, mode);
		return re;
	}
	
	default double[] getCountsForCondition(int condition, ReadCountMode mode) {
		double[] re = new double[getDistinctSequences()];
		for (int i=0; i<re.length; i++)
			re[i] = getCount(i, condition, mode);
		return re;
	}
	
	default int[] getCountsForConditionInt(int condition, ReadCountMode mode) {
		int[] re = new int[getDistinctSequences()];
		for (int i=0; i<re.length; i++)
			re[i] = getCountInt(i, condition, mode);
		return re;
	}
	
	default int[] getCountsForConditionFloor(int condition, ReadCountMode mode) {
		int[] re = new int[getDistinctSequences()];
		for (int i=0; i<re.length; i++)
			re[i] = getCountFloor(i, condition, mode);
		return re;
	}
	
	default NumericArray getCountsForCondition(NumericArray re, int condition, ReadCountMode mode) {
		if (re==null)
			re = NumericArray.createMemory(getDistinctSequences(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getDistinctSequences() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");
		
		for (int distinct=0; distinct<re.length(); distinct++)
			mode.getCount(re,condition,getCount(distinct, condition),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}
	default NumericArray addCountsForCondition(NumericArray re, int condition, ReadCountMode mode) {
		if (re==null)
			re = NumericArray.createMemory(getDistinctSequences(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getDistinctSequences() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");
		
		for (int distinct=0; distinct<re.length(); distinct++)
			mode.addCount(re,condition,getCount(distinct, condition),getMultiplicity(distinct),getWeight(distinct));
		return re;
	}
	
	/**
	 * Gets the counts summed over all conditions for each distinct sequence
	 * @param mode
	 * @return
	 */
	default double[] getTotalCountsForDistincts(ReadCountMode mode) {
		double[] re = new double[getDistinctSequences()];
		for (int d=0; d<re.length; d++)
			for (int i=0; i<getNumConditions(); i++)
				re[d] += getCount(d, i, mode);
		return re;
	}
	
	/**
	 * Gets the counts summed over all distinct sequences for each condition
	 * @param mode
	 * @return
	 */
	default double[] getTotalCountsForConditions(ReadCountMode mode) {
		double[] re = new double[getNumConditions()];
		for (int i=0; i<getDistinctSequences(); i++)
			addCountsForDistinct(i, re, mode);
//			for (int c=0; c<getNumConditions(); c++)
//				re[c] += getCount(i, c, mode);
		return re;
	}
	
	default double getTotalCountOverall(ReadCountMode mode) {
		int re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getTotalCountForDistinct(i, mode);
		return re;
	}
	
	
	
	default int getCountInt(int distinct, int condition, ReadCountMode mode) {
		return mode.computeCountInt(getCount(distinct, condition), getMultiplicity(distinct), getWeight(distinct));
	}
	
	default double getTotalCountForDistinct(int distinct, ReadCountMode mode) {
		if (hasNonzeroInformation()) {
			double re = 0;
			int n = getNonzeroCountIndicesForDistinct(distinct).length;
			for (int i=0; i<n; i++)
				re += getNonzeroCountValue(distinct, i, mode);
			return re;
		}
		
		double re = 0;
		for (int i=0; i<getNumConditions(); i++)
			re += getCount(distinct, i, mode);
		return re;
	}
	
	default double getTotalCountForDistinct(int distinct, ReadCountMode mode, IntDoubleToDoubleFunction op) {
		if (hasNonzeroInformation()) {
			double re = 0;
			int[] inds = getNonzeroCountIndicesForDistinct(distinct);
			for (int i=0; i<inds.length; i++)
				re += op.applyAsDouble(inds[i],getNonzeroCountValue(distinct, i, mode));
			return re;
		}
		
		double re = 0;
		for (int i=0; i<getNumConditions(); i++)
			re += op.applyAsDouble(i,getCount(distinct, i, mode));
		return re;
	}
	
	default int getTotalCountForDistinctInt(int distinct, ReadCountMode mode) {
		if (hasNonzeroInformation()) {
			int re = 0;
			int n = getNonzeroCountIndicesForDistinct(distinct).length;
			for (int i=0; i<n; i++)
				re += getNonzeroCountValueInt(distinct, i, mode);
			return re;
		}
		
		int re = 0;
		for (int i=0; i<getNumConditions(); i++)
			re += getCountInt(distinct, i, mode);
		return re;
	}
	
	default int getTotalCountForConditionInt(int condition, ReadCountMode mode) {
		int re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getCountInt(i, condition, mode);
		return re;
	}
	
	
	/**
	 * Gets the counts summed over all conditions for each distinct sequence
	 * @param mode
	 * @return
	 */
	default int[] getTotalCountsForDistinctsInt(ReadCountMode mode) {
		int[] re = new int[getDistinctSequences()];
		for (int d=0; d<re.length; d++)
			re[d] = getTotalCountForDistinctInt(d, mode);
		return re;
	}
	
	
	/**
	 * Gets the counts summed over all distinct sequences for each condition
	 * @param mode
	 * @return
	 */
	default int[] getTotalCountsForConditionsInt(ReadCountMode mode) {
		int[] re = new int[getNumConditions()];
		for (int c=0; c<getNumConditions(); c++)
			for (int i=0; i<getDistinctSequences(); i++)
				re[c] += getCountInt(i,c,  mode);
		return re;
	}
	
	
	default int getTotalCountOverallInt(ReadCountMode mode) {
		int re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getTotalCountForDistinctInt(i, mode);
		return re;
	}
	
	
	
	default int getCountFloor(int distinct, int condition, ReadCountMode mode) {
		return mode.computeCountFloor(getCount(distinct, condition), getMultiplicity(distinct), getWeight(distinct));
	}
	
	default int getTotalCountForDistinctFloor(int distinct, ReadCountMode mode) {
		if (hasNonzeroInformation()) {
			int re = 0;
			int n = getNonzeroCountIndicesForDistinct(distinct).length;
			for (int i=0; i<n; i++)
				re += getNonzeroCountValueFloor(distinct, i, mode);
			return re;
		}
		
		int re = 0;
		for (int i=0; i<getNumConditions(); i++)
			re += getCountFloor(distinct, i, mode);
		return re;
	}
	
	default int getTotalCountForConditionFloor(int condition, ReadCountMode mode) {
		int re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getCountFloor(i, condition, mode);
		return re;
	}
	
	/**
	 * Gets the counts summed over all conditions for each distinct sequence
	 * @param mode
	 * @return
	 */
	default int[] getTotalCountsForDistinctsFloor(ReadCountMode mode) {
		int[] re = new int[getDistinctSequences()];
		for (int d=0; d<re.length; d++)
			re[d] = getTotalCountForDistinctFloor(d, mode);
		return re;
	}
	
	
	/**
	 * Gets the counts summed over all distinct sequences for each condition
	 * @param mode
	 * @return
	 */
	default int[] getTotalCountsForConditionsFloor(ReadCountMode mode) {
		int[] re = new int[getNumConditions()];
		for (int c=0; c<getNumConditions(); c++)
			for (int i=0; i<getDistinctSequences(); i++)
				re[c] += getCountFloor(i, c, mode);
		return re;
	}
	
	
	default int getTotalCountOverallFloor(ReadCountMode mode) {
		int re = 0;
		for (int i=0; i<getDistinctSequences(); i++)
			re += getTotalCountForDistinctFloor(i, mode);
		return re;
	}
	
	/**
	 * Gets the counts summed over all conditions for each distinct sequence
	 * 
	 * Counts are added to re; re can be null; if re.length or its type do not match, a new one is created
	 * @param re
	 * @param mode
	 * @return
	 */
	default NumericArray getTotalCountsForDistincts(NumericArray re, ReadCountMode mode) {
		if (re!=null)
			re.clear();
		return addTotalCountsForDistincts(re, mode);
	}
	
	default NumericArray addTotalCountsForDistincts(NumericArray re, ReadCountMode mode) {
		if (re==null)
			re = NumericArray.createMemory(getDistinctSequences(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getDistinctSequences() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");

		
		
		for (int d=0; d<getDistinctSequences(); d++)
			if (hasNonzeroInformation()) {
				int[] inds = getNonzeroCountIndicesForDistinct(d);
				for (int i=0; i<inds.length; i++)
					mode.addCount(re,d,getNonzeroCountValueForDistinct(d, i),getMultiplicity(d),getWeight(d));
			}
			else
				for (int c=0; c<getNumConditions(); c++)
					mode.addCount(re,d,getCount(d, c),getMultiplicity(d),getWeight(d));
		return re;
	}
	
	
	/**
	 * 
	 * Gets the counts summed over all distinct sequences for each condition.
	 * 
	 * Counts are added to re; re can be null; if re.length or its type do not match, a new one is created
	 * @param re
	 * @param mode
	 * @return
	 */
	default NumericArray getTotalCountsForConditions(NumericArray re, ReadCountMode mode) {
		if (re!=null)
			re.clear();
		return addTotalCountsForConditions(re, mode);
	}
	
	default NumericArray addTotalCountsForConditions(NumericArray re, ReadCountMode mode) {
		if (re==null)
			re = NumericArray.createMemory(getNumConditions(), NumericArrayType.fromType(mode.getType()));
		else if (re.length()!=getNumConditions() || re.getType().getType()!=mode.getType())
			throw new RuntimeException("Type or length does not match!");

		for (int d=0; d<getDistinctSequences(); d++)
			if (hasNonzeroInformation()){
				int[] inds = getNonzeroCountIndicesForDistinct(d);
				for (int i=0; i<inds.length; i++)
					mode.addCount(re,inds[i],getNonzeroCountValueForDistinct(d, i),getMultiplicity(d),getWeight(d));
			} 
			else
				for (int c=0; c<getNumConditions(); c++)
					mode.addCount(re,c,getCount(d, c),getMultiplicity(d),getWeight(d));
		return re;
	}
	
//	
//	
//	
//	default int getTotalCountUniqueMappings(int condition) {
//		int re = 0;
//		for (int i=0; i<getDistinctSequences(); i++)
//			if (getMultiplicity(i)==1)
//				re+=getCount(i, condition);
//		return re;
//	}
//	default double getTotalCountDivide(int condition) {
//		double re = 0;
//		for (int i=0; i<getDistinctSequences(); i++)
//			re+=(double)getCount(i, condition)/getMultiplicity(i);
//		return re;
//	}
//	default int getTotalCount(int condition) {
//		int re = 0;
//		for (int i=0; i<getDistinctSequences(); i++)
//			re+=getCount(i, condition);
//		return re;
//	}
//	
//
//	
//	default int[] getTotalCount(int[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new int[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] = getTotalCount(i);
//		return re;
//	}
//	
//	default double[] getTotalCountDouble(double[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new double[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] = getTotalCount(i);
//		return re;
//	}
//	
//	default int[] getCounts(int distinct, int[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new int[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] = getCount(distinct,i);
//		return re;
//	}
//	
//	default int[] getCounts(int distinct) {
//		return getCounts(distinct,null);
//	}
//	default int[] getTotalCount() {
//		return getTotalCount(null);
//	}
//	default int[] addTotalCount(int[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new int[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] += getTotalCount(i);
//		return re;
//	}
//	
//	default double[] addTotalCount(double[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new double[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] += getTotalCount(i);
//		return re;
//	}
//	
//	default int[] addCount(int distinct, int[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new int[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] += getCount(distinct, i);
//		return re;
//	}
//	
//	default double[] addCount(int distinct, double[] re, boolean divideByMulti) {
//		if (re==null||re.length!=getNumConditions()) re = new double[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] += (double)getCount(distinct, i)/(divideByMulti?getMultiplicity(distinct):1);
//		return re;
//	}
//	
//	default NumericArray addCount(int distinct, NumericArray re, boolean divideByMulti) {
//		if (re==null||re.length()!=getNumConditions()) re = NumericArray.createMemory(getNumConditions(), divideByMulti?NumericArrayType.Double:NumericArrayType.Integer);
//		for (int i=0; i<getNumConditions(); i++)
//			if (divideByMulti)
//				re.add(i, (double)getCount(distinct, i)/Math.max(1, getMultiplicity(distinct)));
//			else
//				re.add(i,getCount(distinct,i));
//		return re;
//	}
//	
//	
//	default NumericArray addSumCount(NumericArray re, boolean divideByMulti) {
//		for (int d=0; d<getDistinctSequences(); d++) 
//			re = addCount(d, re, divideByMulti);
//		return re;
//	}
//	
//	default NumericArray toNumericArray(boolean divideByMulti) {
//		return addSumCount(null, divideByMulti);
//	}
//	
//	default NumericArray addTotalCount(NumericArray re) {
//		if (re==null||re.length()!=getNumConditions()) re = NumericArray.createMemory(getNumConditions(),NumericArrayType.Integer);
//		for (int i=0; i<getNumConditions(); i++)
//			re.add(i,getTotalCount(i));
//		return re;
//	}
//	
//	default int getSumTotalCount() {
//		int re = 0;
//		for (int i=0; i<getNumConditions(); i++)
//			re += getTotalCount(i);
//		return re;
//	}
//	default int getSumTotalCountUniqueMappings() {
//		int re = 0;
//		for (int i=0; i<getNumConditions(); i++)
//			re += getTotalCountUniqueMappings(i);
//		return re;
//	}
//	default double getSumTotalCountDivide() {
//		double re = 0;
//		for (int i=0; i<getNumConditions(); i++)
//			re += getTotalCountDivide(i);
//		return re;
//	}
//	default int getSumCount(int distinct) {
//		int re = 0;
//		for (int i=0; i<getNumConditions(); i++)
//			re += getCount(distinct,i);
//		return re;
//	}
//	default double getSumCount(int distinct, boolean divideByMulti) {
//		double re = 0;
//		for (int i=0; i<getNumConditions(); i++)
//			re += (double)getCount(distinct,i)/(divideByMulti?getMultiplicity(distinct):1);
//		return re;
//	}
//
//	default int[] getTotalCountForDistinctSequence(int distinct, int[] re) {
//		if (re==null||re.length!=getNumConditions()) re = new int[getNumConditions()];
//		for (int i=0; i<getNumConditions(); i++)
//			re[i] = getCount(distinct,i);
//		return re;
//	}
//	default int[] getTotalCountForDistinctSequence(int distinct) {
//		return getTotalCountForDistinctSequence(distinct,null);
//	}

	/**
	 * Checks whether mismatches are consistent with the genomic sequence and inserts N if not.
	 * @param distinct
	 * @param genomic
	 * @return
	 */
	default CharSequence genomeToRead(int distinct, CharSequence genomic) {
		int c = getVariationCount(distinct);
		if (c==0) return genomic;
		int p = 0;
		CharSequence clip3p = "";
		StringBuilder sb = new StringBuilder();
		
		for (int v=0; v<c; v++) {
			if (isInsertion(distinct, v)) {
				int d = getInsertionPos(distinct, v);
				sb.append(genomic.subSequence(p, d));
				sb.append(getInsertion(distinct, v));
				p = d;
				//throw new RuntimeException("Not tested!");
			} else if (isDeletion(distinct, v)) {
				int d = getDeletionPos(distinct, v);
				sb.append(genomic.subSequence(p, d));
				p = d+getDeletion(distinct, v).length();
				//throw new RuntimeException("Not tested!");
			} else if (isMismatch(distinct, v)) {
				int d = getMismatchPos(distinct, v);
				sb.append(genomic.subSequence(p, d));
				if (getMismatchGenomic(distinct, v).charAt(0)==genomic.charAt(d))
					sb.append(getMismatchRead(distinct, v));
				else
					sb.append('N');
				p = d+1;
			} else if (isSoftclip5p(distinct, v)) {
				sb.append(getSoftclip(distinct, v));
				//throw new RuntimeException("Not tested!");
			} else if (isSoftclip(distinct,v)) {
				clip3p = getSoftclip(distinct, v);
				//throw new RuntimeException("Not tested!");
			}
		}
		sb.append(genomic.subSequence(p,genomic.length()));
		sb.append(clip3p);
		
		return sb.toString();
	}
	
	
	default AlignedReadsVariation getVariation(int distinct, int index) {
		if (isMismatch(distinct, index)) return new AlignedReadsMismatch(getMismatchPos(distinct, index), getMismatchGenomic(distinct, index), getMismatchRead(distinct, index),isVariationFromSecondRead(distinct, index));
		if (isDeletion(distinct, index)) return new AlignedReadsDeletion(getDeletionPos(distinct, index), getDeletion(distinct, index),isVariationFromSecondRead(distinct, index));
		if (isInsertion(distinct, index)) return new AlignedReadsInsertion(getInsertionPos(distinct, index), getInsertion(distinct, index),isVariationFromSecondRead(distinct, index));
		if (isSoftclip(distinct, index)) return new AlignedReadsSoftclip(isSoftclip5p(distinct, index), getSoftclip(distinct, index),isVariationFromSecondRead(distinct, index));
		throw new RuntimeException("Unknown variation type!");
	}
	
	default VarIndel getVarIndel(int distinct, int index) {
		if (isMismatch(distinct, index)) return AlignedReadsDataFactory.createMismatch(getMismatchPos(distinct, index), getMismatchGenomic(distinct, index).charAt(0), getMismatchRead(distinct, index).charAt(0),isVariationFromSecondRead(distinct, index));
		if (isDeletion(distinct, index)) return AlignedReadsDataFactory.createDeletion(getDeletionPos(distinct, index), getDeletion(distinct, index),isVariationFromSecondRead(distinct, index));
		if (isInsertion(distinct, index)) return AlignedReadsDataFactory.createInsertion(getInsertionPos(distinct, index), getInsertion(distinct, index),isVariationFromSecondRead(distinct, index));
		if (isSoftclip(distinct, index)) return AlignedReadsDataFactory.createSoftclip(isSoftclip5p(distinct, index), getSoftclip(distinct, index),isVariationFromSecondRead(distinct, index));
		throw new RuntimeException("Unknown variation type!");
	}
	
	default AlignedReadsVariation[] getVariations(int distinct) {
		AlignedReadsVariation[] re = new AlignedReadsVariation[getVariationCount(distinct)];
		for (int i=0; i<re.length; i++)
			re[i] = getVariation(distinct, i);
		return re;
	}
	
	default ExtendedIterator<VarIndel> getVarIndels(int distinct) {
		return EI.seq(0, getVariationCount(distinct)).map(i->getVarIndel(distinct, i));
	}

	default HasSubreads getHasSubread() {
		if (this instanceof HasSubreads)
			return (HasSubreads) this;
		else if (this instanceof AlignedReadsDataDecorator)
			return ((AlignedReadsDataDecorator) this).getParent().getHasSubread();
		return null;
	}

	default String toString2() {
		StringBuilder sb = new StringBuilder();
		for (int d=0; d<getDistinctSequences(); d++) {
			if (hasId()) 
				sb.append(getId(d)).append(": ");
			if (hasNonzeroInformation()) {
				sb.append("[");
				int[] inds = getNonzeroCountIndicesForDistinct(d);
				for (int i=0; i<inds.length; i++)
					sb.append(inds[i]+":"+getNonzeroCountValueForDistinct(d, i)+",");
				sb.deleteCharAt(sb.length()-1);
				sb.append("]");
			} else
				sb.append(Arrays.toString(getCountsForDistinctInt(d, ReadCountMode.All)));
			sb.append(" x");
			sb.append(getMultiplicity(d));
			if (hasWeights()) 
				sb.append(" (w=").append(String.format("%.2f", getWeight(d))).append(")");
			if (hasGeometry()) 
				sb.append(String.format(" %d|%d|%d", getGeometryBeforeOverlap(d),getGeometryOverlap(d),getGeometryAfterOverlap(d)));
			HasSubreads sr = getHasSubread();
			if (sr!=null) {
				sb.append(" ");
				sr.addSubreadToString(sb, d);
			}
			for (AlignedReadsVariation var : getVariations(d))
				sb.append("\t"+var);
			if (d<getDistinctSequences()-1) sb.append(" ~ ");
		}
		return sb.toString();
	}

	default int hashCode2() {
		int result = 1;
		for (int i=0; i<getDistinctSequences(); i++) {
			int cs = 1;
			for (int j=0; j<getNumConditions(); j++) {
				int e = getCount(i, j);
				
				int elementHash = (int)(e ^ (e>>> 32));
				cs = 31 * cs + elementHash;
			}
			result += cs;
			if (hasId())
				result += Integer.hashCode(getId(i))<<13;
			
			int vs = 1;
			for (int j=0; j<getVariationCount(i); j++) {
				int e = getVariation(i, j).hashCode();
				vs = e;
			}
			result+=31*vs;
		}
		return result;
	}
	
	default boolean equals(Object obj, boolean considerMultiplicity, boolean considerId) {
		return equals(obj,considerMultiplicity,considerId,0);
	}
	default boolean equals(Object obj, boolean considerMultiplicity, boolean considerId, double countTolerance) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (!(obj instanceof AlignedReadsData))
			return false;
		AlignedReadsData other = (AlignedReadsData) obj;
		
		if (getNumConditions()!=other.getNumConditions()) return false;
		if (getDistinctSequences()!=other.getDistinctSequences()) return false;
		
		AlignedReadsVariation[][] a = new AlignedReadsVariation[getDistinctSequences()][];
		AlignedReadsVariation[][] b = new AlignedReadsVariation[getDistinctSequences()][];
		for (int i=0; i<a.length; i++) {
			a[i] = getVariations(i);
			b[i] = other.getVariations(i);
			Arrays.sort(a[i]);
			Arrays.sort(b[i]);
		}
		Ranking<AlignedReadsVariation[]> ranka = new Ranking<AlignedReadsVariation[]>(a,FunctorUtils.arrayComparator()).sort(true);
		Ranking<AlignedReadsVariation[]> rankb = new Ranking<AlignedReadsVariation[]>(b,FunctorUtils.arrayComparator()).sort(true);
		for (int i=0; i<a.length; i++) {
			if (getVariationCount(ranka.getOriginalIndex(i))!=other.getVariationCount(rankb.getOriginalIndex(i))) return false;
			if (considerMultiplicity && getMultiplicity(ranka.getOriginalIndex(i))!=other.getMultiplicity(rankb.getOriginalIndex(i))) return false;
			if (considerId && getId(ranka.getOriginalIndex(i))!=other.getId(rankb.getOriginalIndex(i))) return false;
			if (countTolerance==0) {
				for (int j=0; j<getNumConditions(); j++)
					if (getCount(ranka.getOriginalIndex(i),j)!=other.getCount(rankb.getOriginalIndex(i),j)) return false;
			} else {
				for (int j=0; j<getNumConditions(); j++){
					int ca = getCount(ranka.getOriginalIndex(i),j);
					int cb = other.getCount(rankb.getOriginalIndex(i),j);
					if (Math.abs(ca-cb)>Math.max(ca,cb)*countTolerance) return false;
				}
			}
			
			if (!Arrays.equals(a[i], b[i])) return false;
		}
		
		return true;
	}

	
//	static final byte INT_ID_MAGIC = 0;
//	static final byte LONG_ID_MAGIC = 1;
//	static final byte STRING_ID_MAGIC = 2;
//	static final byte NO_ID_MAGIC = 2;
	
	default void serialize(BinaryWriter out) throws IOException {
		int d = getDistinctSequences();
		int c = getNumConditions();
		
		out.putCInt(d);// distinct
		
		DynamicObject gi = out.getContext().getGlobalInfo();
		if (!gi.hasProperty(CONDITIONSATTRIBUTE))
			out.putCInt(c);//conditions
		
		if (!gi.hasProperty(SPARSEATTRIBUTE) || gi.getEntry(SPARSEATTRIBUTE).asInt()==0) {
			for (int i=0; i<d; i++)
				for (int j=0; j<c; j++)
					out.putCInt(getCount(i, j));
		}
		else if (gi.getEntry(SPARSEATTRIBUTE).asInt()==1){
			int co = 0;
			for (int i=0; i<d; i++)
				for (int j=0; j<c; j++) 
					if (getCount(i, j)>0)
						co++;
			out.putCInt(co);
			int ind = 0;
			for (int i=0; i<d; i++)
				for (int j=0; j<c; j++) {
					if (getCount(i, j)>0) {
						out.putCInt(ind);
						out.putCInt(getCount(i, j));
					}
					ind++;
				}
		}
		else if (gi.getEntry(SPARSEATTRIBUTE).asInt()==2){ // supersparse mode
			for (int i=0; i<d; i++) {
				
				
				if (hasNonzeroInformation()) {
					int[] inds = getNonzeroCountIndicesForDistinct(i);
					out.putCInt(inds.length);
					for (int ii=0; ii<inds.length; ii++) {
						out.putCInt(inds[ii]);
						out.putCInt(getNonzeroCountValueForDistinct(i, ii));
					}
				} else {
					int co = 0;
					for (int j=0; j<c; j++) 
						if (getCount(i, j)>0)
							co++;
					out.putCInt(co);
					for (int j=0; j<c; j++) {
						if (getCount(i,j)>0) {
							out.putCInt(j);
							out.putCInt(getCount(i, j));
						}
					}
				}
			}
		}
		
		for (int i=0; i<d; i++) {
			int v = getVariationCount(i);
			out.putCInt(v);
			for (int j=0; j<v; j++) {
				CharSequence ch;
				if (isMismatch(i, j)) {
					out.putCShort(DefaultAlignedReadsData.encodeMismatch(getMismatchPos(i, j), getMismatchGenomic(i, j).charAt(0), getMismatchRead(i, j).charAt(0),isVariationFromSecondRead(i, j)));
					ch = DefaultAlignedReadsData.encodeMismatchIndel(getMismatchPos(i, j), getMismatchGenomic(i, j).charAt(0), getMismatchRead(i, j).charAt(0));
				} else if (isDeletion(i, j)) {
					out.putCShort(DefaultAlignedReadsData.encodeDeletion(getDeletionPos(i, j), getDeletion(i, j),isVariationFromSecondRead(i, j)));
					ch = DefaultAlignedReadsData.encodeDeletionIndel(getDeletionPos(i, j), getDeletion(i, j));
					out.putCInt(ch.length());
				} else if (isInsertion(i, j)) {
					out.putCShort(DefaultAlignedReadsData.encodeInsertion(getInsertionPos(i, j), getInsertion(i, j),isVariationFromSecondRead(i, j)));
					ch = DefaultAlignedReadsData.encodeInsertionIndel(getInsertionPos(i, j), getInsertion(i, j));
					out.putCInt(ch.length());
				} else if (isSoftclip(i, j)) {
					out.putCShort(DefaultAlignedReadsData.encodeSoftclip(isSoftclip5p(i, j), getSoftclip(i, j),isVariationFromSecondRead(i, j)));
					ch = DefaultAlignedReadsData.encodeSoftclipSequence(isSoftclip5p(i, j), getSoftclip(i, j));
					out.putCInt(ch.length());
				} else 
					throw new RuntimeException("Neither mismatch nor deletion nor insertion!");
				out.putAsciiChars(ch);
			}
		}
			
		for (int i=0; i<d; i++)
			out.putCInt(getMultiplicity(i));
		
		
		if (out.getContext().getGlobalInfo().getEntry(AlignedReadsData.HASIDATTRIBUTE).asInt()==1) {
			for (int i=0; i<d; i++)
				out.putCInt(getId(i));
		} 
		
		if (out.getContext().getGlobalInfo().getEntry(AlignedReadsData.HASWEIGHTATTRIBUTE).asInt()==1) {
			for (int i=0; i<d; i++)
				out.putFloat(getWeight(i));
		} 
		
		if (out.getContext().getGlobalInfo().getEntry(AlignedReadsData.HASGEOMETRYATTRIBUTE).asInt()==1) {
			for (int i=0; i<d; i++)
				out.putCInt(DefaultAlignedReadsData.encodeGeometry(getGeometryBeforeOverlap(i),getGeometryOverlap(i),getGeometryAfterOverlap(i)));
		} 
		
//		if (hasIntId()) {
//			out.put(INT_ID_MAGIC);
//			for (int i=0; i<d; i++)
//				out.putCInt(getIntId(i));
//		}
//		else if (hasLongId()) {
//			out.put(LONG_ID_MAGIC);
//			for (int i=0; i<d; i++)
//				out.putCLong(getLongId(i));
//		}
//		else if (hasId()) {
//			out.put(STRING_ID_MAGIC);
//			for (int i=0; i<d; i++)
//				out.putString(getId(i));
//		}
//		else
//			out.put(NO_ID_MAGIC);

		if (this instanceof HasBarcodes) {
			HasBarcodes bc = (HasBarcodes)this;
			if (!gi.hasProperty(BARCODEATTRIBUTE)) {
				out.putCInt(bc.getBarcodeLength());			
			}
			
			for (int i=0; i<d; i++) {
				
				if (hasNonzeroInformation())
					for (int j : getNonzeroCountIndicesForDistinct(i)) {
						DnaSequence[] bcs = bc.getBarcodes(i, j);
						for (int b=0; b<bcs.length;b++)
							FileUtils.writeBitVector(bcs[b], out, false);
					}
				else
					for (int j=0; j<c; j++) {
						if (getCount(i,j)>0) {
							DnaSequence[] bcs = bc.getBarcodes(i, j);
							for (int b=0; b<bcs.length;b++)
								FileUtils.writeBitVector(bcs[b], out, false);
						}
					}
				
				
			}
		}
		
		
		if (this instanceof HasSubreads) {
			HasSubreads bc = (HasSubreads)this;
			
			for (int i=0; i<d; i++) {
				out.putCInt(bc.getNumSubreads(i));
				for (int s=0; s<bc.getNumSubreads(i); s++) {
					out.putCInt(bc.getSubreadId(i, s));
					if (s<bc.getNumSubreads(i)-1) 
						out.putCInt(bc.getSubreadEnd(i, s,-1));
				}
				out.putCInt(bc.getGapPositions(i).countInt());
				IntIterator it = bc.getGapPositions(i);
				while (it.hasNext()) out.putCInt(it.nextInt());
			}
		}
		
	}
	default int getNumConditionsWithCounts() {
		int re = 0;
		for (int i=0; i<getNumConditions(); i++)
			if (getTotalCountForCondition(i, ReadCountMode.All)>0)
				re++;
		return re;
	}
	
	
	/**
	 * The default implementation produces a subread that either has {@link SingleEndSubreadsSemantic} or {@link PairedEndSubreadsSemantic}.
	 * 
	 * Single end:
	 * If this subread is supposed to be a sense read (parameter=true), the ids will be as follows:
	 * ###########
	 * 00000000000
	 * If this subread is supposed to be an antisense read (parameter=false), the ids will be as follows:
	 * ##########
	 * 1111111111
	 * 
	 * 
	 * Paired end:
	 * If this subread is supposed to be a sense read (parameter=true), the ids will be as follows:
	 * ##########
	 *      ###########
	 * 0000011111222222
	 * 
	 * ##########
	 *           ###########
	 * 000000000022222222222
	 * ##########
	 *       ####
	 * 0000001111
	 *
	 * etc
	 * 
	 * If this subread is supposed to be an antisense read (parameter=false), the ids will be as follows:
	 * ##########
	 *      ###########
	 * 2222211111000000
	 * 
	 * ##########
	 *           ###########
	 * 222222222200000000000
	 *
	 * ##########
	 *       ####
	 * 2222221111
	 * 
	 * etc
	 * @param sense
	 * @return
	 */
	default HasSubreads asSubreads(boolean sense) {
		return new HasSubreads() {
			@Override
			public int getNumSubreads(int distinct) {
				if (!hasGeometry()) return 1;
				int re = 0;
				if (getGeometryBeforeOverlap(distinct)>0) re++;
				if (getGeometryOverlap(distinct)>0) re++;
				if (getGeometryAfterOverlap(distinct)>0) re++;
				return re;
			}

			@Override
			public int getSubreadStart(int distinct, int index) {
				if (index==0) return 0;
				
				if (sense && index==1) {
					if (getGeometryBeforeOverlap(distinct)>0) return getGeometryBeforeOverlap(distinct);
					return getGeometryOverlap(distinct);
				}
				if (index==1) {
					if (getGeometryAfterOverlap(distinct)>0) return getGeometryAfterOverlap(distinct);
					return getGeometryOverlap(distinct);
				}
				
				if (sense) 
					return getGeometryBeforeOverlap(distinct)+getGeometryOverlap(distinct);
				
				return getGeometryAfterOverlap(distinct)+getGeometryOverlap(distinct);

			}

			@Override
			public int getSubreadId(int distinct, int index) {
				if (!hasGeometry()) return sense?0:1;
				
				if (sense && index==0) {
					if (getGeometryBeforeOverlap(distinct)>0) return 0;
					if (getGeometryOverlap(distinct)>0) return 1;
					return 2;
				}
				if (sense && index==1) {
					if (getGeometryBeforeOverlap(distinct)>0 && getGeometryOverlap(distinct)>0) return 1;
					return 2;
				}
				if (sense) {
					return 2;
				}
				if (index==0) {
					if (getGeometryAfterOverlap(distinct)>0) return 0;
					if (getGeometryOverlap(distinct)>0) return 1;
					return 2;
				}
				if (index==1) {
					if (getGeometryAfterOverlap(distinct)>0 && getGeometryOverlap(distinct)>0) return 1;
					return 2;
				}
				return 2;
			}

			@Override
			public IntIterator getGapPositions(int distinct) {
				if (!hasGeometry() || getGeometryOverlap(distinct)>0)
					return IntIterator.empty;
				return IntIterator.singleton(getGeometryBeforeOverlap(distinct));
			}
		};
	}

	default boolean isAnyAmbigousMapping() {
		for (int i=0; i<getDistinctSequences(); i++)
			if (isAmbigousMapping(i)) return true;
		return false;
	}
	
	default boolean isAnyUniqueMapping() {
		for (int i=0; i<getDistinctSequences(); i++)
			if (!isAmbigousMapping(i)) return true;
		return false;
	}

	
	
	default ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> iterateDistinct(ImmutableReferenceGenomicRegion<?> rgr) {
		return EI.seq(0, getDistinctSequences()).map(d->new ImmutableReferenceGenomicRegion<>(rgr.getReference(),rgr.getRegion(),new SelectDistinctSequenceAlignedReadsData(this, d)));
	}

	default ExtendedIterator<AlignedReadsData> iterateDistinct() {
		return EI.seq(0, getDistinctSequences()).map(d->new SelectDistinctSequenceAlignedReadsData(this, d));
	}

	default ImmutableReferenceGenomicRegion<SelectDistinctSequenceAlignedReadsData> restrictToDistinct(ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, int... distinct) {
		return new ImmutableReferenceGenomicRegion<>(read.getReference(),read.getRegion(),restrictToDistinct(distinct));
	}

	default SelectDistinctSequenceAlignedReadsData restrictToDistinct(int... distinct) {
		return new SelectDistinctSequenceAlignedReadsData(this, distinct);
	}

}

