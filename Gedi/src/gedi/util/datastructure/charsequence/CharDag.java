package gedi.util.datastructure.charsequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
import java.util.function.BiConsumer;
import java.util.function.IntFunction;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.StringUtils;
import gedi.util.datastructure.tree.Trie;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

public class CharDag {

	
	
	private CharIterator text;
	
	private SequenceVariant[] variants;
	
	
	public CharDag(CharIterator text, ExtendedIterator<SequenceVariant> variants) {
		this.text = text;
		this.variants = variants.toArray(SequenceVariant.class);
		Arrays.sort(this.variants);
	}
	
	public SequenceVariant[] getVariants() {
		return variants;
	}
	
	
	public <R> void traverse(CharDagVisitor<R> v, BiConsumer<R,Iterator<SequenceVariant>> re) {
		traverse(i->i==0?v:null, re);
	}
	
	public <R> void traverse(IntFunction<CharDagVisitor<R>> v, BiConsumer<R,Iterator<SequenceVariant>> re) {
		
		int nextVariantIndex = 0;
		int pos = -1;
		ArrayList<CharDagVisitor<R>> vs = new ArrayList<CharDagVisitor<R>>();
		TreeMap<Integer, ArrayList<CharDagVisitor<R>>> cont = new TreeMap<>();
		
		while (text.hasNext()) {
			pos++;
			
			CharDagVisitor<R> vvv = v.apply(pos);
			if (vvv!=null)
				vs.add(vvv);
			
			while (nextVariantIndex<variants.length && variants[nextVariantIndex].pos==pos) {
				
				ArrayList<CharDagVisitor<R>> bvs = new ArrayList<CharDagVisitor<R>>();
				for (CharDagVisitor<R> vv : vs) 
					bvs.add(vv.branch());
			
				if (bvs.size()>1) 
					Collections.sort(bvs);

//				for (CharDagVisitor<R> sv : bvs) 
//					sv.shift(pos,-variants[nextVariantIndex].from.length, variants[nextVariantIndex]);
				
				for (int pp=0; pp<variants[nextVariantIndex].to.length; pp++) {
					char c = variants[nextVariantIndex].to[pp];
					
					int index = 0;
//					bvs.get(0).shift(pos,1, variants[nextVariantIndex]);
					Iterator<R> r = bvs.get(0).accept(c,variants[nextVariantIndex],pp);
					accept(re,r, bvs.get(0),variants[nextVariantIndex]);
					for (int i=1; i<bvs.size(); i++) {
//						bvs.get(i).shift(pos,1, variants[nextVariantIndex]);
						r = bvs.get(i).accept(c,pos);
						if (bvs.get(i).compareTo(bvs.get(index))!=0) {
							index++;
							accept(re,r, bvs.get(i),variants[nextVariantIndex]);					
						}
						bvs.set(index, bvs.get(i));
					}
					while (bvs.size()!=index+1) {
						bvs.remove(bvs.size()-1);
					}
				}
				for (CharDagVisitor<R> sv : bvs) 
					sv.addVariant(variants[nextVariantIndex]);
				
				int conti = pos+variants[nextVariantIndex].from.length;
				ArrayList<CharDagVisitor<R>> cl = cont.get(conti);
				if (cl==null) cont.put(conti, bvs);
				else cl.addAll(bvs);
				nextVariantIndex++;
			}
			
			ArrayList<CharDagVisitor<R>> cl = cont.remove(pos);
			if (cl!=null) 
				vs.addAll(cl);
			
			if (vs.size()>1) 
				Collections.sort(vs);
			
			char c = text.nextChar();
			
			int index = 0;
			Iterator<R> r = vs.get(0).accept(c,pos);
			accept(re, r, vs.get(0));					
			for (int i=1; i<vs.size(); i++) {
				r = vs.get(i).accept(c,pos);
				if (vs.get(i).compareTo(vs.get(index))!=0) {
					index++;
					accept(re, r, vs.get(i));
				}
				vs.set(index, vs.get(i));
			}
			while (vs.size()!=index+1)
				vs.remove(vs.size()-1);
		
			for (CharDagVisitor<R> vv : vs)
				vv.prune(pos);
			
		}
		
	}
	
	
	
	private <R> void accept(BiConsumer<R,Iterator<SequenceVariant>> re, Iterator<R> r, CharDagVisitor<R> v) {
		if (r!=null) {
			while (r.hasNext())
				re.accept(r.next(),v.getVariants());					
		}
	}
	private <R> void accept(BiConsumer<R,Iterator<SequenceVariant>> re, Iterator<R> r, CharDagVisitor<R> v, SequenceVariant currentVar) {
		if (r!=null) {
			while (r.hasNext())
				re.accept(r.next(),EI.wrap(v.getVariants()).chain(EI.singleton(currentVar)));					
		}
	}


	public static interface CharDagVisitor<R> extends Comparable<CharDagVisitor<R>> {
		Iterator<R> accept(char c, int pos);
		Iterator<R> accept(char c, SequenceVariant variant, int varToPos);
		Iterator<SequenceVariant> getVariants();
		CharDagVisitor<R> branch();
//		void shift(int pos, int shift, SequenceVariant variant);
		void addVariant(SequenceVariant variant);
		void prune(int pos);
	}
	
	
	public static class KmerCharDagVisitor implements CharDagVisitor<CharSequence> {

		private CharRingBuffer buff;
		private int n = 0;
		private LinkedList<SequenceVariant> correction = new LinkedList<SequenceVariant>();
		
		private KmerCharDagVisitor() {
		}
			
		public KmerCharDagVisitor(int k) {
			this.buff = new CharRingBuffer(k);
		}

		@Override
		public Iterator<CharSequence> accept(char c, int pos) {
			this.buff.add(c);
			if (++n>=this.buff.capacity())
				return EI.singleton(this.buff.toString()+"@"+pos+", "+StringUtils.toString(correction));
			return null;
		}
		
		@Override
		public Iterator<CharSequence> accept(char c, SequenceVariant var, int varToPos) {
			this.buff.add(c);
			if (++n>=this.buff.capacity())
				return EI.singleton(this.buff.toString()+"@"+varToPos+" in "+var+", "+StringUtils.toString(correction));
			return null;
		}
		
		@Override
		public int compareTo(CharDagVisitor<CharSequence> o) {
			return buff.compareTo(((KmerCharDagVisitor)o).buff);
		}

		@Override
		public CharDagVisitor<CharSequence> branch() {
			KmerCharDagVisitor re = new KmerCharDagVisitor();
			re.buff = new CharRingBuffer(buff);
			re.n = n;
			re.correction.addAll(correction);
			return re;
		}

//		@Override
//		public void shift(int pos, int shift, SequenceVariant variant) {
//			if (correction.size()>0 && correction.getLast().getPos()==pos && correction.getLast().getVariant().equals(variant)) {
//				correction.getLast().incrementLen(shift);
//				if (correction.getLast().getLen()==0) correction.removeLast();
//			} else
//				correction.add(new PosLengthVariant(pos, shift,variant));
//		}
		
		@Override
		public void addVariant(SequenceVariant variant) {
			correction.add(variant);
		}

		@Override
		public void prune(int pos) {
			Iterator<SequenceVariant> it = correction.iterator();
			while (it.hasNext() && it.next().getEndPositionWithFrom()<pos-buff.capacity())
				it.remove();			
		}

		@Override
		public Iterator<SequenceVariant> getVariants() {
			return EI.wrap(correction);
		}
		
	}
	
	public static class SequenceVariant implements Comparable<SequenceVariant> {
		int pos;
		char[] from;
		char[] to;
		String label;
		
		public SequenceVariant(int pos, char[] from, char[] to, String label) {
			this.pos = pos;
			this.from = from;
			this.to = to;
			this.label = label;
		}
		
		public int getPosition() {
			return pos;
		}
		
		public int getEndPositionWithFrom() {
			return pos+from.length;
		}
		
		public int getLengthDelta() {
			return to.length-from.length;
		}
		
		public int getToLength() {
			return to.length;
		}
		
		public int getFromLength() {
			return from.length;
		}

		@Override
		public int compareTo(SequenceVariant o) {
			return Integer.compare(pos, o.pos);
		}

		public int end() {
			return pos+from.length;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + Arrays.hashCode(from);
			result = prime * result + pos;
			result = prime * result + Arrays.hashCode(to);
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
			SequenceVariant other = (SequenceVariant) obj;
			if (!Arrays.equals(from, other.from))
				return false;
			if (pos != other.pos)
				return false;
			if (!Arrays.equals(to, other.to))
				return false;
			return true;
		}

		@Override
		public String toString() {
			if (label!=null)
				return label;
			return pos+" "+String.valueOf(from)+">"+String.valueOf(to);
		}

		public SequenceVariant compact() {
			int pref;
			for (pref=0; pref<to.length && pref<from.length && to[pref]==from[pref]; pref++);
			int suff;
			for (suff=0; suff<to.length-pref && suff<from.length-pref && to[to.length-1-suff]==from[from.length-1-suff]; suff++);
			
			if (pref>0 || suff>0) {
				pos+=pref;
				from = Arrays.copyOfRange(from, pref, from.length-suff);
				to = Arrays.copyOfRange(to, pref, to.length-suff);
			}
			return this;
		}
		
		
		
	}
	
	
	public static void main(String[] args) {
//		//            0123456789
//		String seq = "TTAAAAATGAGA";
//		//            TTAAGAGA 
//		Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = new Trie<List<ImmutableReferenceGenomicRegion<String>>>();
//		aho.put("K", new ArrayList<>());
//		aho.put(SequenceUtils.translate("TTA"), new ArrayList<>());
//
//		CharDag dag = new CharDag(CharIterator.fromCharSequence(seq),EI.wrap(
//				new SequenceVariant(4, "AAAT".toCharArray(), "G".toCharArray(), "AAAG->G")
//				));
////		,EI.wrap(
////				new SequenceVariant(2, "2".toCharArray(), "2a".toCharArray(), "3->a3"),
////				new SequenceVariant(3, "3".toCharArray(), "a3".toCharArray(), "2->2a")
////				));
//		
//		dag.traverse(i->i<3?new TranslateCharDagVisitor<>(aho.ahoCorasickVisitor(l->l*3),SequenceUtils.codeTrie,false):null,(res, varit)->{
//			System.out.println(res);
//			System.out.println(EI.wrap(varit).concat(" "));
//		});
		
		
		String seq = "0123456789";
		Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = new Trie<List<ImmutableReferenceGenomicRegion<String>>>();
		aho.put("0a", new ArrayList<>());

		CharDag dag = new CharDag(CharIterator.fromCharSequence(seq),EI.wrap(
				new SequenceVariant(1, "123".toCharArray(), "a".toCharArray(), "mut")
				));
//		,EI.wrap(
//				new SequenceVariant(2, "2".toCharArray(), "2a".toCharArray(), "3->a3"),
//				new SequenceVariant(3, "3".toCharArray(), "a3".toCharArray(), "2->2a")
//				));
		
		dag.traverse(aho.ahoCorasickVisitor(),(res, varit)->{
			System.out.println(res);
			System.out.println(EI.wrap(varit).concat(" "));
		});
	}
	
	
	
}
