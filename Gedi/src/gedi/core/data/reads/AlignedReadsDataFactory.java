package gedi.core.data.reads;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.function.UnaryOperator;

import cern.colt.bitvector.BitVector;
import gedi.util.ArrayUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.collections.doublecollections.DoubleArrayList;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.mutable.MutablePair;
import gedi.util.sequence.DnaSequence;

public class AlignedReadsDataFactory {

	protected int currentDistinct;
	protected int conditions;
	
	protected ArrayList<int[]> count = new ArrayList<int[]>();
	protected ArrayList<int[]> nonzeros;
	protected IntArrayList nextcount;
	
//	protected ArrayList<IntArrayList> var = new ArrayList<IntArrayList>();
//	protected ArrayList<ArrayList<CharSequence>> indels = new ArrayList<ArrayList<CharSequence>>();
	protected ArrayList<TreeSet<VarIndel>> var = new ArrayList<TreeSet<VarIndel>>();
	protected ArrayList<MutablePair<int[],int[]>> subs = new ArrayList<MutablePair<int[],int[]>>();
	
	protected IntArrayList multiplicity = new IntArrayList();
	protected DoubleArrayList weights = new DoubleArrayList();
	protected IntArrayList geometry = new IntArrayList();
	protected IntArrayList ids = new IntArrayList();
//	protected ArrayList ids = new ArrayList();
	
	protected boolean hasSubs = false;
	
	protected ArrayList<ArrayList<DnaSequence>[]> barcodes = new ArrayList<ArrayList<DnaSequence>[]>();
	protected boolean hasBarcodes = false;
	
	public AlignedReadsDataFactory(int conditions) {
		this(conditions,conditions>5);
	}
	public AlignedReadsDataFactory(int conditions, boolean sparse) {
		this.conditions = conditions;
		if (sparse) {
			nonzeros = new ArrayList<int[]>();
			nextcount = new IntArrayList();
		}
	}
	
	
	
	public static void removeId(DefaultAlignedReadsData ard) {
		ard.ids = null;
	}
	
	public AlignedReadsDataFactory start() {
		hasSubs = false;
		currentDistinct=-1;
		count.clear();
		if (isSparse()) {
			nonzeros.clear();
			nextcount.clear();
		}
		var.clear();
		subs.clear();
		multiplicity.clear();
		weights.clear();
		geometry.clear();
		ids.clear();
		barcodes.clear();
		hasBarcodes = false;
		return this;
	}
	public int getDistinctSequences() {
		return currentDistinct+1;
	}
	public AlignedReadsDataFactory newDistinctSequence() {
		currentDistinct++;
		setMultiplicity(1);
		
		this.count.add(new int[conditions]);
		if (isSparse()) {
			this.nonzeros.add(new int[conditions]);
			this.nextcount.add(0);
		}
		this.var.add(new TreeSet<VarIndel>());
		this.subs.add(new MutablePair<int[],int[]>());
		this.barcodes.add(new ArrayList[conditions]);
		
		return this;
	}
	private void checkDistinct() {
		if (currentDistinct<0) throw new RuntimeException("Call newDistinctSequence first!");
	}
	private void checkDistinct(int distinct) {
		if (currentDistinct<0 || distinct>currentDistinct) throw new RuntimeException("Call newDistinctSequence first!");
	}
	
//	private void checkPos(int pos) {
//		if (var.get(currentDistinct).size()==0) return;
//		int lastPos = DefaultAlignedReadsData.pos(var.get(currentDistinct).getLastInt());
//		if (lastPos>=pos) throw new RuntimeException("Provide variations in-order: "+lastPos+" >= "+pos);
//	}

	public AlignedReadsDataFactory setMultiplicity(int m) {
		checkDistinct();
		multiplicity.set(currentDistinct, m);
		return this;
	}
	
	public AlignedReadsDataFactory setId(int id) {
		checkDistinct();
		if (ids.size()<currentDistinct)
			throw new RuntimeException("Call setId for all distinct sequences!");
		setId(currentDistinct, id);
		return this;
	}	
	public AlignedReadsDataFactory setWeight(float weight) {
		checkDistinct();
		if (weights.size()<currentDistinct)
			throw new RuntimeException("Call setWeight for all distinct sequences!");
		setWeight(currentDistinct, weight);
		return this;
	}
	
	public AlignedReadsDataFactory setGeometry(int before, int overlap, int after) {
		checkDistinct();
		if (geometry.size()<currentDistinct)
			throw new RuntimeException("Call setGeometry for all distinct sequences!");
		setGeometry(currentDistinct, before, overlap, after);
		return this;
	}
	
	public AlignedReadsDataFactory setCount(int[] count) {
		return setCount(currentDistinct,count);
	}
	public AlignedReadsDataFactory setCount(int distinct,int[] count) {
		checkDistinct(distinct);
		if (isSparse()) {
			for (int i=0; i<count.length; i++)
				if (count[i]>0)
					setCount(distinct,i,count[i]);
		} 
		else 
			System.arraycopy(count, 0, this.count.get(distinct), 0, conditions);
		return this;
	}
	
	public AlignedReadsDataFactory setCount(int condition, int count) {
		return setCount(currentDistinct, condition, count, null);
	}
	public AlignedReadsDataFactory setCount(int distinct,int condition, int count) {
		return setCount(distinct, condition, count, null);
	}
		
	public AlignedReadsDataFactory setCount(int condition, int count, DnaSequence[] barcodes) {
		return setCount(currentDistinct, condition, count,barcodes);
	}
	public AlignedReadsDataFactory setCount(int distinct,int condition, int count, DnaSequence[] barcodes) {
		if (count<=0) throw new RuntimeException("Count is non-positive!");
		if (barcodes!=null && barcodes.length!=count) 
			throw new RuntimeException("Barcodes and counts are inconsistent!");
		if (barcodes!=null)
			hasBarcodes = true;
		
		checkDistinct(distinct);
		
		
		if (isSparse()) {
			int nc = nextcount.getInt(distinct);
			int[] nz = nonzeros.get(distinct);
			int[] co = this.count.get(distinct);
			if (nc>0 && nz[nc-1]>condition)
				throw new RuntimeException("Internal error: Violated contract (must insert by increasing conditions!)!");
			
			boolean update = nc>0 && nz[nc-1]==condition;
			
			if (update) {
				co[nc-1] =  count;
				if (barcodes!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					bc[nc-1].clear();
					bc[nc-1].addAll(Arrays.asList(barcodes));
				}
			} else {
				co[nc] =  count;
				nz[nc] =  condition;
				if (barcodes!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					bc[nc] = new ArrayList<DnaSequence>();
					bc[nc].addAll(Arrays.asList(barcodes));
				}
				this.nextcount.increment(distinct);
			}
			
		} else {
			this.count.get(distinct)[condition]=count;
			if (barcodes!=null)
				this.barcodes.get(distinct)[condition]=new ArrayList<DnaSequence>(Arrays.asList(barcodes));
		}
		return this;
	}
	

	public AlignedReadsDataFactory incrementCountSingle(int condition) {
		return incrementCountSingle(currentDistinct, condition);
	}
	public AlignedReadsDataFactory incrementCountSingle(int distinct,int condition) {
		return incrementCountSingle(distinct,condition,null);
	}
	public AlignedReadsDataFactory incrementCountSingle(int distinct,int condition, DnaSequence barcode) {
		if (barcode!=null)
			hasBarcodes = true;
		
		checkDistinct(distinct);
		
		if (isSparse()) {
			int nc = nextcount.getInt(distinct);
			int[] nz = nonzeros.get(distinct);
			int[] co = this.count.get(distinct);
			if (nc>0 && nz[nc-1]>condition)
				throw new RuntimeException("Internal error: Violated contract (must insert by increasing conditions!)!");
			
			boolean update = nc>0 && nz[nc-1]==condition;
			
			if (update) {
				co[nc-1]++;
				if (barcode!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					if (bc[nc-1]==null) throw new RuntimeException("Either add barcodes always or never!");
					bc[nc-1].add(barcode);
				}
			} else {
				co[nc]++;
				nz[nc] =  condition;
				if (barcode!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					if (bc[nc]!=null) throw new RuntimeException("Either add barcodes always or never!");
					bc[nc] = new ArrayList<DnaSequence>();
					bc[nc].add(barcode);
				}
				this.nextcount.increment(distinct);
			}
			
		} else {
			this.count.get(distinct)[condition]++;
			if (barcode!=null) {
				if (this.barcodes.get(distinct)[condition]==null)
					this.barcodes.get(distinct)[condition]=new ArrayList<DnaSequence>();
				this.barcodes.get(distinct)[condition].add(barcode);
			}
		}
		return this;
	}

	
	public AlignedReadsDataFactory incrementCount(int condition, int count) {
		return incrementCount(currentDistinct, condition, count);
	}
	public AlignedReadsDataFactory incrementCount(int distinct,int condition, int count) {
		return incrementCount(distinct,condition,count,null);
	}
	public AlignedReadsDataFactory incrementCount(int distinct,int condition, int count, DnaSequence[] barcodes) {
		if (count<0) throw new RuntimeException("Count is negative!");
		if (barcodes!=null && barcodes.length!=count) 
			throw new RuntimeException("Barcodes and counts are inconsistent!");
		if (barcodes!=null)
			hasBarcodes = true;
		if (count==0) return this;
		
		checkDistinct(distinct);
		
		if (isSparse()) {
			int nc = nextcount.getInt(distinct);
			int[] nz = nonzeros.get(distinct);
			int[] co = this.count.get(distinct);
			if (nc>0 && nz[nc-1]>condition)
				throw new RuntimeException("Internal error: Violated contract (must insert by increasing conditions!)!");
			
			boolean update = nc>0 && nz[nc-1]==condition;
			
			if (update) {
				co[nc-1] +=  count;
				if (barcodes!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					if (bc[nc-1]==null) throw new RuntimeException("Either add barcodes always or never!");
					bc[nc-1].addAll(Arrays.asList(barcodes));
				}
			} else {
				co[nc] +=  count;
				nz[nc] =  condition;
				if (barcodes!=null) {
					ArrayList<DnaSequence>[] bc = this.barcodes.get(distinct);
					if (bc[nc]!=null) throw new RuntimeException("Either add barcodes always or never!");
					bc[nc] = new ArrayList<DnaSequence>();
					bc[nc].addAll(Arrays.asList(barcodes));
				}
				this.nextcount.increment(distinct);
			}
			
		} else {
			this.count.get(distinct)[condition]+=count;
			if (barcodes!=null) {
				if (this.barcodes.get(distinct)[condition]==null)
					this.barcodes.get(distinct)[condition]=new ArrayList<DnaSequence>();
				this.barcodes.get(distinct)[condition].addAll(Arrays.asList(barcodes));
			}
		}
		return this;
	}
	
	
	public AlignedReadsDataFactory incrementCount(NumericArray a) {
		checkDistinct();
		for (int i=0; i<getNumConditions(); i++)
			incrementCount(i,a.getInt(i));
		return this;
	}
	
//	public AlignedReadsDataFactory addBarcode(int condition, DnaSequence[] barcodes) {
//		checkDistinct();
//		if (isSparse()) {
//			throw new RuntimeException("Barcodes are not stored in a sparse array! Fix it and dont forget makeDistinct!");
//		}			
//		if (this.barcodes.get(currentDistinct)[condition]==null)
//			this.barcodes.get(currentDistinct)[condition]=new ArrayList<>();
//		this.barcodes.get(currentDistinct)[condition].addAll(Arrays.asList(barcodes));
//		hasBarcodes=true;
//		return this;
//	}
	
//	public AlignedReadsDataFactory addBarcode(int condition, String barcode) {
//		checkDistinct();
//		if (isSparse()) {
//			throw new RuntimeException("Barcodes are not stored in a sparse array! Fix it and dont forget makeDistinct! If they are saved to a CIT, they might be sparsely stored (if there are more than 5 conditions), and that's fine (and good) for further procedures)");
//		}			
//		if (this.barcodes.get(currentDistinct)[condition]==null)
//			this.barcodes.get(currentDistinct)[condition]=new ArrayList<>();
//		this.barcodes.get(currentDistinct)[condition].add(new DnaSequence(barcode));
//		hasBarcodes=true;
//		return this;
//	}
	
//	public AlignedReadsDataFactory setBarcodes(ArrayList<DnaSequence>[] barcodes) {
//		checkDistinct();
//		this.barcodes.set(currentDistinct,barcodes);
//		hasBarcodes=true;
//		return this;
//	}

	public AlignedReadsDataFactory setId(int distinct,int id) {
		checkDistinct(distinct);
		ids.set(distinct, id);
		return this;
	}
	public AlignedReadsDataFactory setWeight(int distinct,float weight) {
		checkDistinct(distinct);
		weights.set(distinct, weight);
		return this;
	}
	public AlignedReadsDataFactory incrementWeight(int distinct,float weight) {
		checkDistinct(distinct);
		weights.set(distinct, weights.getDouble(distinct)+weight);
		return this;
	}
	
	public AlignedReadsDataFactory setGeometry(int distinct,int before, int overlap, int after) {
		checkDistinct(distinct);
		geometry.set(distinct, DefaultAlignedReadsData.encodeGeometry(before, overlap, after));
		return this;
	}
	
//	
//	public AlignedReadsDataFactory setId(int distinct,long id) {
//		checkDistinct(distinct);
//		checkId(Long.class);
//		ids.set(distinct, id);
//		return this;
//	}
//	
//	public AlignedReadsDataFactory setId(int distinct,String id) {
//		checkDistinct(distinct);
//		checkId(String.class);
//		ids.set(distinct, id);
//		return this;
//	}
//	
//	
//	private void checkId(Class<?> cls) {
//		if (currentDistinct>0) {
//			if (ids.get(0).getClass()!=cls)
//				throw new RuntimeException("Id types must be homogeneous!");
//		}
//	}

	public AlignedReadsDataFactory setMultiplicity(int distinct,int m) {
		checkDistinct(distinct);
		multiplicity.set(distinct, m);
		return this;
	}
	
	public int getMultiplicity(int distinct) {
		checkDistinct(distinct);
		return multiplicity.getInt(distinct);
	}
	
	
//	public AlignedReadsDataFactory addBarcode(int distinct, int condition, String barcode) {
//		checkDistinct(distinct);
//		if (isSparse()) {
//			throw new RuntimeException("Barcodes are not stored in a sparse array! Fix it!");
//		}		
//		
//		if (this.barcodes.get(distinct)[condition]==null)
//			this.barcodes.get(distinct)[condition]=new ArrayList<>();
//		this.barcodes.get(distinct)[condition].add(new DnaSequence(barcode));
//		hasBarcodes=true;
//		return this;
//	}

	public AlignedReadsDataFactory setSubread(int[] sub, int[] gap) {
		hasSubs = true;
		checkDistinct();
		MutablePair<int[], int[]> p = this.subs.get(currentDistinct);
		p.Item1 = sub;
		p.Item2 = gap;
		return this;
	}
	
	public AlignedReadsDataFactory addMismatch(int pos, char genomic, char read, boolean isFromSecondRead) {
		getCurrentVariationBuffer().add(createMismatch(pos, genomic, read,isFromSecondRead));
		return this;
	}
	
	public AlignedReadsDataFactory addSoftclip(boolean p5, CharSequence read, boolean isFromSecondRead) {
		if (read.length()>0)
			getCurrentVariationBuffer().add(createSoftclip(p5, read,isFromSecondRead));
		return this;
	}
	
	public AlignedReadsDataFactory addDeletion(int pos, CharSequence genomic, boolean isFromSecondRead) {
		getCurrentVariationBuffer().add(createDeletion(pos, genomic,isFromSecondRead));
		return this;
	}
	
	public AlignedReadsDataFactory addInsertion(int pos, CharSequence read, boolean isFromSecondRead) {
		getCurrentVariationBuffer().add(createInsertion(pos, read,isFromSecondRead));
		return this;
	}
	
	protected Collection<VarIndel> getCurrentVariationBuffer() {
		checkDistinct();
		return this.var.get(currentDistinct);
	}
	
	public VarIndel createVarIndel(String variationString) {
		
		String rest = variationString.substring(1);
		int n = StringUtils.countPrefixInt(rest);
		if (n==0) throw new IllegalArgumentException("No pos found in "+variationString);
		
		int pos = Integer.parseInt(rest.substring(0, n));
		rest = rest.substring(n);
		
		boolean second = rest.endsWith("r");
		if (second)
			rest = rest.substring(0,rest.length()-1);
		
		if (variationString.startsWith("M"))
			return createMismatch(pos, rest.charAt(0), rest.charAt(1),second);
		
		if (variationString.startsWith("I"))
			return createInsertion(pos, rest,second);
		
		if (variationString.startsWith("D"))
			return createDeletion(pos, rest,second);
		
		throw new IllegalArgumentException("Must start with M/I/D: "+variationString);
	}
	
	public void addVariation(AlignedReadsVariation vari) {
		if (vari.isDeletion())
			addDeletion(vari.getPosition(), vari.getReferenceSequence(),vari.isFromSecondRead());
		else if (vari.isInsertion())
			addInsertion(vari.getPosition(), vari.getReadSequence(),vari.isFromSecondRead());
		else if (vari.isMismatch())
			addMismatch(vari.getPosition(), vari.getReferenceSequence().charAt(0),vari.getReadSequence().charAt(0),vari.isFromSecondRead());
		else if (vari.isSoftclip())
			addSoftclip(vari.getPosition()==0, vari.getReadSequence(),vari.isFromSecondRead());
	}

	public void addVariationToOtherRead(AlignedReadsVariation vari) {
		if (vari.isDeletion())
			addDeletion(vari.getPosition(), SequenceUtils.getDnaComplement(vari.getReferenceSequence()),!vari.isFromSecondRead());
		else if (vari.isInsertion())
			addInsertion(vari.getPosition(), SequenceUtils.getDnaComplement(vari.getReadSequence()),!vari.isFromSecondRead());
		else if (vari.isMismatch())
			addMismatch(vari.getPosition(), SequenceUtils.getDnaComplement(vari.getReferenceSequence().charAt(0)),SequenceUtils.getDnaComplement(vari.getReadSequence().charAt(0)),!vari.isFromSecondRead());
		else if (vari.isSoftclip())
			addSoftclip(vari.getPosition()==0, SequenceUtils.getDnaComplement(vari.getReadSequence()),!vari.isFromSecondRead());
	}

	public static VarIndel createVariation(AlignedReadsVariation vari) {
		if (vari.isDeletion())
			return createDeletion(vari.getPosition(), vari.getReferenceSequence(),vari.isFromSecondRead());
		else if (vari.isInsertion())
			return createInsertion(vari.getPosition(), vari.getReadSequence(),vari.isFromSecondRead());
		else if (vari.isMismatch())
			return createMismatch(vari.getPosition(), vari.getReferenceSequence().charAt(0),vari.getReadSequence().charAt(0),vari.isFromSecondRead());
		else if (vari.isSoftclip())
			return createSoftclip(vari.getPosition()==0, vari.getReadSequence(),vari.isFromSecondRead());
		throw new RuntimeException("Unknown variation type!");
	}
	
	public static VarIndel createSoftclip(boolean p5, CharSequence read, boolean onSecondRead) {
		VarIndel v = new VarIndel();
		v.var = DefaultAlignedReadsData.encodeSoftclip(p5, read,onSecondRead);
		v.indel = DefaultAlignedReadsData.encodeSoftclipSequence(p5, read);
		return v;
	}
	
	
	public static VarIndel createMismatch(int pos, char genomic, char read, boolean onSecondRead) {
		VarIndel v = new VarIndel();
		v.var = DefaultAlignedReadsData.encodeMismatch(pos, genomic, read,onSecondRead);
		v.indel = DefaultAlignedReadsData.encodeMismatchIndel(pos, genomic, read);
		return v;
	}
	
	public static VarIndel createDeletion(int pos, CharSequence genomic, boolean onSecondRead) {
		VarIndel v = new VarIndel();
		v.var = DefaultAlignedReadsData.encodeDeletion(pos, genomic,onSecondRead);
		v.indel = DefaultAlignedReadsData.encodeDeletionIndel(pos, genomic);
		return v;
	}
	
	public static VarIndel createInsertion(int pos, CharSequence read, boolean onSecondRead) {
		VarIndel v = new VarIndel();
		v.var = DefaultAlignedReadsData.encodeInsertion(pos, read,onSecondRead);
		v.indel = DefaultAlignedReadsData.encodeInsertionIndel(pos, read);
		return v;
	}
	
	public static class VarIndel implements AlignedReadsVariation {
		public short var;
		public CharSequence indel;
		@Override
		public int compareTo(AlignedReadsVariation o) {
			VarIndel a = (VarIndel) o;
			int re = Integer.compare(DefaultAlignedReadsData.pos(var), DefaultAlignedReadsData.pos(a.var));
			if (re==0) re = Short.compare(var, a.var);
			if (re==0) re = StringUtils.compare(indel, a.indel);
			return re;
		}
		
		public int getPosition() {
			return DefaultAlignedReadsData.pos(var);
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((indel == null) ? 0 : indel.hashCode());
			result = prime * result + var;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			VarIndel other = (VarIndel) obj;
			if (indel == null) {
				if (other.indel != null)
					return false;
			} else if (!indel.equals(other.indel))
				return false;
			if (var != other.var)
				return false;
			return true;
		}
		
		@Override
		public int getType() {
			return DefaultAlignedReadsData.type(var);
		}

		
		@Override
		public String toString() {
			if (DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_MISMATCH) 
				return new AlignedReadsMismatch(DefaultAlignedReadsData.pos(var), indel.subSequence(0, 1), indel.subSequence(1, 2),DefaultAlignedReadsData.isSecondRead(var)).toString();
			if (DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_DELETION) 
				return new AlignedReadsDeletion(DefaultAlignedReadsData.pos(var), indel,DefaultAlignedReadsData.isSecondRead(var)).toString();
			if (DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_INSERTION) 
				return new AlignedReadsInsertion(DefaultAlignedReadsData.pos(var), indel,DefaultAlignedReadsData.isSecondRead(var)).toString();
			if (DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_SOFTCLIP) 
				return new AlignedReadsSoftclip(DefaultAlignedReadsData.pos(var)==0, indel,DefaultAlignedReadsData.isSecondRead(var)).toString();

			return "VarIndel [var=" + var + ", indel=" + indel + "]";
		}


		@Override
		public boolean isFromSecondRead() {
			return DefaultAlignedReadsData.isSecondRead(var);
		}


		@Override
		public boolean isMismatch() {
			return DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_MISMATCH;
		}


		@Override
		public boolean isDeletion() {
			return DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_DELETION;
		}


		@Override
		public boolean isInsertion() {
			return DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_INSERTION;
		}


		@Override
		public boolean isSoftclip() {
			return DefaultAlignedReadsData.type(var)==DefaultAlignedReadsData.TYPE_SOFTCLIP;
		}


		@Override
		public CharSequence getReferenceSequence() {
			if (isSoftclip() || isInsertion()) return "";
			if (isMismatch()) indel.subSequence(0, 1);
			return indel;
		}


		@Override
		public CharSequence getReadSequence() {
			if (isSoftclip() || isDeletion()) return "";
			if (isMismatch()) indel.subSequence(1, 2);
			return indel;
		}


		@Override
		public VarIndel reposition(int newPos) {
			if (newPos>DefaultAlignedReadsData.MAX_POSITION) return null;
			var = DefaultAlignedReadsData.encodePos(newPos,var);
			return this;
		}

		public VarIndel complement() {
			indel = SequenceUtils.getDnaComplement(indel);
			return this;
		}
		
	}
	
	private static class GeomDistincter {
		TreeSet<VarIndel> vars;
		int geometry;
		public GeomDistincter(TreeSet<VarIndel> vars, int geometry) {
			super();
			this.vars = vars;
			this.geometry = geometry;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + geometry;
			result = prime * result + ((vars == null) ? 0 : vars.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			GeomDistincter other = (GeomDistincter) obj;
			if (geometry != other.geometry)
				return false;
			if (vars == null) {
				if (other.vars != null)
					return false;
			} else if (!vars.equals(other.vars))
				return false;
			return true;
		}

		
	}
	
	private static class SubDistincter {
		TreeSet<VarIndel> vars;
		MutablePair<int[],int[]> subs;
		public SubDistincter(TreeSet<VarIndel> vars, MutablePair<int[], int[]> subs) {
			super();
			this.vars = vars;
			this.subs = subs;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((subs == null) ? 0 : Arrays.hashCode(subs.Item1));
			result = prime * result + ((subs == null) ? 0 : Arrays.hashCode(subs.Item2));
			result = prime * result + ((vars == null) ? 0 : vars.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			SubDistincter other = (SubDistincter) obj;
			if (subs == null) {
				if (other.subs != null)
					return false;
			} else if (!Arrays.equals(subs.Item1, other.subs.Item1) || !Arrays.equals(subs.Item2, other.subs.Item2))
				return false;
			if (vars == null) {
				if (other.vars != null)
					return false;
			} else if (!vars.equals(other.vars))
				return false;
			return true;
		}
		
		
	}
	
	public boolean areTrulyDistinct() {
		if (hasSubs) return areTrulyDistinctSubs();
		
		HashSet<GeomDistincter> h = new HashSet<GeomDistincter>();
		for (int i=0; i<var.size(); i++)
			h.add(new GeomDistincter(var.get(i), geometry.getInt(i)));
		return h.size()==var.size();
		
//		HashSet<TreeSet<VarIndel>> h = new HashSet<TreeSet<VarIndel>>();
//		for (TreeSet<VarIndel> t : var)
//			h.add(t);
//		return h.size()==var.size();
	}
	private boolean areTrulyDistinctSubs() {
		
		HashSet<SubDistincter> h = new HashSet<SubDistincter>();
		for (int i=0; i<var.size(); i++)
			h.add(new SubDistincter(var.get(i), subs.get(i)));
		return h.size()==var.size();
	}
	
	public void makeDistinct() {
		if (hasSubs) {
			makeDistinctSubs();
			return;
		}
		makeDistinctGeom();
		
//		HashMap<TreeSet<VarIndel>,Integer> h = new HashMap<TreeSet<VarIndel>,Integer>();
//		for (int i=0; i<var.size(); i++) {
//			Integer idx = h.get(var.get(i));
//			if (idx==null) h.put(var.get(i), i);
//			else {
//				if (!isSparse()) {
//					ArrayUtils.add(count.get(idx), count.get(i));
//					for (int c=0; c<conditions; c++) {
//						if (barcodes.get(idx)[c]==null && barcodes.get(i)[c]!=null)
//							barcodes.get(idx)[c] = new ArrayList<DnaSequence>();
//						if (barcodes.get(i)[c]!=null)
//							barcodes.get(idx)[c].addAll(barcodes.get(i)[c]);
//					}
//				} else {
//					IntArrayList cc = new IntArrayList();
//					ArrayList<ArrayList<DnaSequence>> bc = hasBarcodes?new ArrayList<ArrayList<DnaSequence>>():null;
//					IntArrayList cn = new IntArrayList();
//					int idxi = 0;
//					int ii = 0;
//					while (idxi<nextcount.getInt(idx) && ii<nextcount.getInt(i)) {
//						if (nonzeros.get(idx)[idxi]==nonzeros.get(i)[ii]) { 
//							cn.add(nonzeros.get(idx)[idxi]);
//							cc.add(count.get(idx)[idxi]+count.get(i)[ii]);
//							if (hasBarcodes) {
//								ArrayList<DnaSequence> l = barcodes.get(idx)[idxi];
//								l.addAll(barcodes.get(i)[ii]);
//								bc.add(l);
//							}
//							ii++;
//							idxi++;
//						}
//						else if (nonzeros.get(idx)[idxi]<nonzeros.get(i)[ii]){
//							cn.add(nonzeros.get(idx)[idxi]);
//							cc.add(count.get(idx)[idxi]);
//							if (hasBarcodes) {
//								bc.add(barcodes.get(idx)[idxi]);
//							}
//							idxi++;
//						}
//						else {
//							cn.add(nonzeros.get(i)[ii]);
//							cc.add(count.get(i)[ii]);
//							if (hasBarcodes) {
//								bc.add(barcodes.get(i)[ii]);
//							}
//							ii++;
//						}
//					}
//					while (idxi<nextcount.getInt(idx)) {
//						cn.add(nonzeros.get(idx)[idxi]);
//						cc.add(count.get(idx)[idxi]);
//						if (hasBarcodes) {
//							bc.add(barcodes.get(idx)[idxi]);
//						}
//						idxi++;
//					}
//					while (ii<nextcount.getInt(i)) {
//						cn.add(nonzeros.get(i)[ii]);
//						cc.add(count.get(i)[ii]);
//						if (hasBarcodes) {
//							bc.add(barcodes.get(i)[ii]);
//						}
//						ii++;
//					}
//					count.set(idx,cc.toIntArray());
//					if (hasBarcodes) barcodes.set(idx, bc.toArray(new ArrayList[0]));
//					nonzeros.set(idx, cn.toIntArray());
//					nextcount.set(idx,cc.size());
//				}
//				if (ids.size()>0)
//					ids.set(idx, Math.min(ids.getInt(idx),ids.getInt(i)));
//			}
//		}
//		
//		
//		
//		BitVector keep = new BitVector(var.size());
//		for (Integer i : h.values()) 
//			keep.putQuick(i, true);
//		
//		count = restrict(count,keep);
//		if (isSparse()) {
//			nonzeros = restrict(nonzeros,keep);
//			nextcount = restrict(nextcount,keep);
//		}
//		var = restrict(var,keep);
//		subs = restrict(subs,keep);
//		multiplicity = restrict(multiplicity,keep);
//		weights = restrict(weights,keep);
//		geometry = restrict(geometry,keep);
//		ids = restrict(ids,keep);
//		barcodes = restrict(barcodes,keep);
//		
//		currentDistinct = count.size()-1;
	}
	
	private void makeDistinctGeom() {
		
		HashMap<GeomDistincter,Integer> h = new HashMap<GeomDistincter,Integer>();
		for (int i=0; i<var.size(); i++) {
			GeomDistincter dist = new GeomDistincter(var.get(i), geometry.getInt(i));
			
			Integer idx = h.get(dist);
			if (idx==null) h.put(dist, i);
			else {
				if (!isSparse()) {
					ArrayUtils.add(count.get(idx), count.get(i));
					for (int c=0; c<conditions; c++) {
						if (barcodes.get(idx)[c]==null && barcodes.get(i)[c]!=null)
							barcodes.get(idx)[c] = new ArrayList<DnaSequence>();
						if (barcodes.get(i)[c]!=null)
							barcodes.get(idx)[c].addAll(barcodes.get(i)[c]);
					}
				} else {
					IntArrayList cc = new IntArrayList();
					ArrayList<ArrayList<DnaSequence>> bc = hasBarcodes?new ArrayList<ArrayList<DnaSequence>>():null;
					IntArrayList cn = new IntArrayList();
					int idxi = 0;
					int ii = 0;
					while (idxi<nextcount.getInt(idx) && ii<nextcount.getInt(i)) {
						if (nonzeros.get(idx)[idxi]==nonzeros.get(i)[ii]) { 
							cn.add(nonzeros.get(idx)[idxi]);
							cc.add(count.get(idx)[idxi]+count.get(i)[ii]);
							if (hasBarcodes) {
								ArrayList<DnaSequence> l = barcodes.get(idx)[idxi];
								l.addAll(barcodes.get(i)[ii]);
								bc.add(l);
							}
							ii++;
							idxi++;
						}
						else if (nonzeros.get(idx)[idxi]<nonzeros.get(i)[ii]){
							cn.add(nonzeros.get(idx)[idxi]);
							cc.add(count.get(idx)[idxi]);
							if (hasBarcodes) {
								bc.add(barcodes.get(idx)[idxi]);
							}
							idxi++;
						}
						else {
							cn.add(nonzeros.get(i)[ii]);
							cc.add(count.get(i)[ii]);
							if (hasBarcodes) {
								bc.add(barcodes.get(i)[ii]);
							}
							ii++;
						}
					}
					while (idxi<nextcount.getInt(idx)) {
						cn.add(nonzeros.get(idx)[idxi]);
						cc.add(count.get(idx)[idxi]);
						if (hasBarcodes) {
							bc.add(barcodes.get(idx)[idxi]);
						}
						idxi++;
					}
					while (ii<nextcount.getInt(i)) {
						cn.add(nonzeros.get(i)[ii]);
						cc.add(count.get(i)[ii]);
						if (hasBarcodes) {
							bc.add(barcodes.get(i)[ii]);
						}
						ii++;
					}
					count.set(idx,cc.toIntArray());
					if (hasBarcodes) barcodes.set(idx, bc.toArray(new ArrayList[0]));
					nonzeros.set(idx, cn.toIntArray());
					nextcount.set(idx,cc.size());
				}
				if (ids.size()>0)
					ids.set(idx, Math.min(ids.getInt(idx),ids.getInt(i)));
			}
		}
		
		
		
		BitVector keep = new BitVector(var.size());
		for (Integer i : h.values()) 
			keep.putQuick(i, true);
		
		count = restrict(count,keep);
		if (isSparse()) {
			nonzeros = restrict(nonzeros,keep);
			nextcount = restrict(nextcount,keep);
		}
		var = restrict(var,keep);
		subs = restrict(subs,keep);
		multiplicity = restrict(multiplicity,keep);
		weights = restrict(weights,keep);
		geometry = restrict(geometry,keep);
		ids = restrict(ids,keep);
		barcodes = restrict(barcodes,keep);
		
		currentDistinct = count.size()-1;
	}

	private void makeDistinctSubs() {
	
	HashMap<SubDistincter,Integer> h = new HashMap<SubDistincter,Integer>();
	for (int i=0; i<var.size(); i++) {
		SubDistincter dist = new SubDistincter(var.get(i), subs.get(i));
		
		Integer idx = h.get(dist);
		if (idx==null) h.put(dist, i);
		else {
			if (!isSparse()) {
				ArrayUtils.add(count.get(idx), count.get(i));
				for (int c=0; c<conditions; c++) {
					if (barcodes.get(idx)[c]==null && barcodes.get(i)[c]!=null)
						barcodes.get(idx)[c] = new ArrayList<DnaSequence>();
					if (barcodes.get(i)[c]!=null)
						barcodes.get(idx)[c].addAll(barcodes.get(i)[c]);
				}
			} else {
				IntArrayList cc = new IntArrayList();
				ArrayList<ArrayList<DnaSequence>> bc = hasBarcodes?new ArrayList<ArrayList<DnaSequence>>():null;
				IntArrayList cn = new IntArrayList();
				int idxi = 0;
				int ii = 0;
				while (idxi<nextcount.getInt(idx) && ii<nextcount.getInt(i)) {
					if (nonzeros.get(idx)[idxi]==nonzeros.get(i)[ii]) { 
						cn.add(nonzeros.get(idx)[idxi]);
						cc.add(count.get(idx)[idxi]+count.get(i)[ii]);
						if (hasBarcodes) {
							ArrayList<DnaSequence> l = barcodes.get(idx)[idxi];
							l.addAll(barcodes.get(i)[ii]);
							bc.add(l);
						}
						ii++;
						idxi++;
					}
					else if (nonzeros.get(idx)[idxi]<nonzeros.get(i)[ii]){
						cn.add(nonzeros.get(idx)[idxi]);
						cc.add(count.get(idx)[idxi]);
						if (hasBarcodes) {
							bc.add(barcodes.get(idx)[idxi]);
						}
						idxi++;
					}
					else {
						cn.add(nonzeros.get(i)[ii]);
						cc.add(count.get(i)[ii]);
						if (hasBarcodes) {
							bc.add(barcodes.get(i)[ii]);
						}
						ii++;
					}
				}
				while (idxi<nextcount.getInt(idx)) {
					cn.add(nonzeros.get(idx)[idxi]);
					cc.add(count.get(idx)[idxi]);
					if (hasBarcodes) {
						bc.add(barcodes.get(idx)[idxi]);
					}
					idxi++;
				}
				while (ii<nextcount.getInt(i)) {
					cn.add(nonzeros.get(i)[ii]);
					cc.add(count.get(i)[ii]);
					if (hasBarcodes) {
						bc.add(barcodes.get(i)[ii]);
					}
					ii++;
				}
				count.set(idx,cc.toIntArray());
				if (hasBarcodes) barcodes.set(idx, bc.toArray(new ArrayList[0]));
				nonzeros.set(idx, cn.toIntArray());
				nextcount.set(idx,cc.size());
			}
			if (ids.size()>0)
				ids.set(idx, Math.min(ids.getInt(idx),ids.getInt(i)));
		}
	}
	
	
	
	BitVector keep = new BitVector(var.size());
	for (Integer i : h.values()) 
		keep.putQuick(i, true);
	
	count = restrict(count,keep);
	if (isSparse()) {
		nonzeros = restrict(nonzeros,keep);
		nextcount = restrict(nextcount,keep);
	}
	var = restrict(var,keep);
	subs = restrict(subs,keep);
	multiplicity = restrict(multiplicity,keep);
	weights = restrict(weights,keep);
	geometry = restrict(geometry,keep);
	ids = restrict(ids,keep);
	barcodes = restrict(barcodes,keep);
	
	currentDistinct = count.size()-1;
}
	
	private <T> ArrayList<T> restrict(ArrayList<T> l, BitVector keep) {
		ArrayList<T> re = new ArrayList<T>();
		for (int i=0; i<l.size(); i++)
			if (keep.getQuick(i))
				re.add(l.get(i));
		return re;
	}
	private IntArrayList restrict(IntArrayList l, BitVector keep) {
		IntArrayList re = new IntArrayList();
		for (int i=0; i<l.size(); i++)
			if (keep.getQuick(i))
				re.add(l.getInt(i));
		return re;
	}
	private DoubleArrayList restrict(DoubleArrayList l, BitVector keep) {
		DoubleArrayList re = new DoubleArrayList();
		for (int i=0; i<l.size(); i++)
			if (keep.getQuick(i))
				re.add(l.getDouble(i));
		return re;
	}
	
	public static DefaultAlignedReadsData createSimple(int[] count,boolean sparse) {
		DefaultAlignedReadsData re = new DefaultAlignedReadsData();
		re.conditions = count.length;
		if (sparse) {
			int nz = 0;
			for (int c : count) if (c>0) nz++;
			re.count = new int[1][nz];
			re.nonzeros = new int[1][nz];
			int i = 0;
			for (int ci=0; ci<count.length; ci++)
				if (count[ci]>0) {
					re.count[0][i] = count[ci];
					re.nonzeros[0][i++] = ci;
				}
		}
		else
			re.count = new int[][] {count};
		re.var = new short[1][0];
		re.indels = new CharSequence[1][0];
		re.multiplicity = new int[] {1};
		return re;
	}

	private DefaultAlignedReadsData fill(DefaultAlignedReadsData re) {
		re.conditions = getNumConditions();
		re.count = convInt(count,nextcount);
		if (isSparse()) 
			re.nonzeros = convInt(nonzeros,nextcount);
		re.var = convShort(var);
		re.indels = convS(var);
		re.multiplicity = multiplicity.toIntArray();
		if (ids.size()==re.count.length)
			re.ids = ids.toIntArray();
		else if (ids.size()>0)
			throw new RuntimeException("Call setId for each distinct sequence or for none!");
		
		if (weights.size()==re.count.length)
			re.weights = weights.toFloatArray();
		else if (weights.size()>0)
			throw new RuntimeException("Call setWeight for each distinct sequence or for none!");
		
		if (geometry.size()==re.count.length)
			re.geometry = geometry.toIntArray();
		else if (geometry.size()>0)
			throw new RuntimeException("Call setGeometry for each distinct sequence or for none!");
		
		return re;
	}

	public DefaultAlignedReadsData createDefaultOrBarcode() {
		if (hasBarcodes)
			return createBarcode();
		else
			return create();
	}
	
	public <C extends AlignedReadsData> C create(Class<C> cls) {
		if (cls.equals(DefaultAlignedReadsData.class))
			return (C) create();
		if (cls.equals(BarcodedAlignedReadsData.class))
			return (C) createBarcode();
		throw new RuntimeException("Cannot automatically create class "+cls);
	}
	
	public DefaultAlignedReadsData create() {
		return fill(new DefaultAlignedReadsData());
	}
	
	public BarcodedAlignedReadsData createBarcode() {
		return createBarcode(new BarcodedAlignedReadsData());
	}
	
	public SubreadsAlignedReadsData createSubread() {
		SubreadsAlignedReadsData re = new SubreadsAlignedReadsData();
		fill(re);
		re.subreadGeom = new int[subs.size()][];
		re.gaps = new int[subs.size()][];
		for (int i=0; i<subs.size(); i++) {
			re.subreadGeom[i] = subs.get(i).Item1;
			re.gaps[i] = subs.get(i).Item2;
		}
		return re;
	}
	public BarcodedAlignedReadsData createBarcode(BarcodedAlignedReadsData re) {
		fill(re);
		
		if (isSparse()) {
			re.barcodes  = new DnaSequence[currentDistinct+1][][];
			for (int d=0; d<re.barcodes.length; d++){
				re.barcodes[d] = new DnaSequence[re.nonzeros[d].length][];
				for (int c=0; c<re.getNonzeroCountIndicesForDistinct(d).length; c++) {
					re.barcodes[d][c] = barcodes.get(d)[c].toArray(new DnaSequence[0]);
					if (re.barcodes[d][c].length!=re.count[d][c])
						throw new RuntimeException("Barcodes and counts do not match!");
				}
			}
		} else {
			re.barcodes  = new DnaSequence[currentDistinct+1][conditions][];
			for (int d=0; d<re.barcodes.length; d++){
				for (int c=0; c<re.barcodes[d].length; c++) {
					re.barcodes[d][c] = barcodes.get(d)[c]==null?new DnaSequence[0]:barcodes.get(d)[c].toArray(new DnaSequence[0]);
					if (re.barcodes[d][c].length!=re.count[d][c])
						throw new RuntimeException("Barcodes and counts do not match!");
				}
			}
		}
		return re;
	}
	
	public AlignedReadsData createDigital() {
		DigitalAlignedReadsData re = new DigitalAlignedReadsData();
		if (isSparse() || count.size()!=1 || var.get(0).size()>0) throw new RuntimeException();
		
		re.count = new BitVector(count.get(0).length);
		for (int i=0; i<re.count.size(); i++)
			re.count.putQuick(i, count.get(0)[i]>0);
		return re;
	}


	private CharSequence[][] convS(ArrayList<TreeSet<VarIndel>> var) {
		CharSequence[][] re = new CharSequence[currentDistinct+1][];
		for (int i=0; i<re.length; i++) {
			re[i] = new CharSequence[var.get(i).size()];
			int j=0;
			for (VarIndel v : var.get(i))
				re[i][j++] = v.indel;
		}
		return re;
	}

	private int[][] convInt(ArrayList<int[]> a) {
		int[][] re = new int[currentDistinct+1][];
		for (int i=0; i<re.length; i++) 
			re[i] = a.get(i).clone();
		return re;
	}
	
	private int[][] convInt(ArrayList<int[]> a, IntArrayList nextInd) {
		int[][] re = new int[currentDistinct+1][];
		for (int i=0; i<re.length; i++) {
			if (nextInd==null)
				re[i] = a.get(i).clone();
			else {
				re[i] = new int[nextInd.getInt(i)];
				System.arraycopy(a.get(i), 0, re[i], 0, nextInd.getInt(i));
			}
		}
		return re;
	}

	private short[][] convShort(ArrayList<TreeSet<VarIndel>> var) {
		short[][] re = new short[currentDistinct+1][];
		for (int i=0; i<re.length; i++) {
			re[i] = new short[var.get(i).size()];
			int j=0;
			for (VarIndel v : var.get(i))
				re[i][j++] = v.var;
		}
		return re;
	}


	public int getNumConditions() {
		return conditions;
	}
	public boolean isSparse() {
		return nonzeros!=null;
	}
	
	
	public AlignedReadsDataFactory add(AlignedReadsData ard) {
		for (int d=0; d<ard.getDistinctSequences(); d++)
			add(ard,d);
		return this;
	}
	
	public AlignedReadsDataFactory add(AlignedReadsData ard, UnaryOperator<AlignedReadsVariation> varPred) {
		for (int d=0; d<ard.getDistinctSequences(); d++)
			add(ard,d,varPred);
		return this;
	}
	
	/**
	 * 
	 * @param ard
	 * @param distinct take the distinct from ard!
	 */
	public AlignedReadsDataFactory add(AlignedReadsData ard, int distinct) {
		newDistinctSequence();
		BarcodedAlignedReadsData bc = ard instanceof BarcodedAlignedReadsData?(BarcodedAlignedReadsData)ard:null;
		
		if (ard.hasNonzeroInformation()) {
			int[] is = ard.getNonzeroCountIndicesForDistinct(distinct);
			for (int ii=0; ii<is.length; ii++) {
				int i=is[ii];
				setCount(i, ard.getNonzeroCountValueForDistinct(distinct, ii),bc==null?null:bc.getNonZeroBarcodes(distinct, ii));
			}
		}
		else {
			for (int i=0; i<ard.getNumConditions(); i++)
				if (ard.getCount(distinct, i)>0)
					setCount(i, ard.getCount(distinct, i),bc==null?null:bc.getBarcodes(distinct, i));
		}
		
		setMultiplicity(ard.getMultiplicity(distinct));
		for (int i=0; i<ard.getVariationCount(distinct); i++)
			addVariation(ard.getVariation(distinct, i));
		if (ard.hasId())
			setId(ard.getId(distinct));
		if (ard.hasWeights())
			setWeight(ard.getWeight(distinct));
		if (ard.hasGeometry())
			setGeometry(ard.getGeometryBeforeOverlap(distinct), ard.getGeometryOverlap(distinct), ard.getGeometryAfterOverlap(distinct));
//		if (ard instanceof BarcodedAlignedReadsData) {
//			BarcodedAlignedReadsData bc = (BarcodedAlignedReadsData)ard;
//			for (int i=0; i<ard.getNumConditions(); i++)
//				for (int b=0; b<ard.getCount(distinct, i); b++)
//					addBarcode(i, bc.getBarcodes(distinct, i)[b].toString());
//		}
		return this;
	}
	
	public AlignedReadsDataFactory add(AlignedReadsData ard, int distinct, UnaryOperator<AlignedReadsVariation> varPred) {
		return add(ard,distinct,varPred,false);
	}
	
	public AlignedReadsDataFactory add(AlignedReadsData ard, int distinct, UnaryOperator<AlignedReadsVariation> varPred, boolean removeGeometry) {
		return add(ard,distinct,varPred,removeGeometry,true);
	}
	/**
	 * if varpred is null, no variation is added!
	 * @param ard
	 * @param distinct take the distinct from ard!
	 */
	public AlignedReadsDataFactory add(AlignedReadsData ard, int distinct, UnaryOperator<AlignedReadsVariation> varPred, boolean removeGeometry, boolean setcounts) {
		newDistinctSequence();
		BarcodedAlignedReadsData bc = ard instanceof BarcodedAlignedReadsData?(BarcodedAlignedReadsData)ard:null;
		if (setcounts) {
			if (ard.hasNonzeroInformation()) {
				int[] is = ard.getNonzeroCountIndicesForDistinct(distinct);
				for (int ii=0; ii<is.length; ii++) {
					int i=is[ii];
					setCount(i, ard.getNonzeroCountValueForDistinct(distinct, ii),bc==null?null:bc.getNonZeroBarcodes(distinct, ii));
				}
			}
			else {
				for (int i=0; i<ard.getNumConditions(); i++)
					if (ard.getCount(distinct, i)>0)
						setCount(i, ard.getCount(distinct, i),bc==null?null:bc.getBarcodes(distinct, i));
			}
		}
		setMultiplicity(ard.getMultiplicity(distinct));
		if (varPred!=null)
			for (int i=0; i<ard.getVariationCount(distinct); i++) {
				AlignedReadsVariation vari = ard.getVariation(distinct, i);
				vari = varPred.apply(vari);
				if (vari!=null)
					addVariation(vari);
			}
		if (ard.hasId())
			setId(ard.getId(distinct));
		if (ard.hasWeights())
			setWeight(ard.getWeight(distinct));
		if (!removeGeometry && ard.hasGeometry())
			setGeometry(ard.getGeometryBeforeOverlap(distinct), ard.getGeometryOverlap(distinct), ard.getGeometryAfterOverlap(distinct));
//		if (ard instanceof BarcodedAlignedReadsData) {
//			BarcodedAlignedReadsData bc = (BarcodedAlignedReadsData)ard;
//			if (ard.hasNonzeroInformation()) {
//				int[] is = ard.getNonzeroCountIndicesForDistinct(distinct);
//				for (int ii=0; ii<is.length; ii++) {
//					int i=is[ii];
//					setBarcodes(i, bc.getNonZeroBarcodes(distinct, ii));
//				}
//			}
//			else {
//				for (int i=0; i<ard.getNumConditions(); i++)
//					if (bc.getBarcodes(distinct, i).length>0)
//						setBarcodes(i, bc.getBarcodes(distinct, i));
//			}
//			
////			for (int i=0; i<ard.getNumConditions(); i++)
////				for (int b=0; b<ard.getCount(distinct, i); b++)
////					addBarcode(i, bc.getBarcodes(distinct, i)[b].toString());
//		}
		if (ard instanceof SubreadsAlignedReadsData) {
			SubreadsAlignedReadsData sc = (SubreadsAlignedReadsData) ard;
			setSubread(sc.subreadGeom[distinct],sc.gaps[distinct]);
		}
		return this;
	}

	public void incrementCount(AlignedReadsData ard) {
		for (int i=0; i<ard.getNumConditions(); i++)
			incrementCount(i, ard.getTotalCountForConditionInt(i,ReadCountMode.All));
	}

	

	
	
	
}
