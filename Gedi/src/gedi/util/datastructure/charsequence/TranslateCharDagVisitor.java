package gedi.util.datastructure.charsequence;

import java.util.Iterator;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.charsequence.CharDag.CharDagVisitor;
import gedi.util.datastructure.charsequence.CharDag.SequenceVariant;
import gedi.util.datastructure.tree.Trie;

public class TranslateCharDagVisitor<R> implements CharDagVisitor<R> {

	private char[] codon = new char[3];
	private int p;
	private CharDagVisitor<R> parent;
	private boolean revert;
	private Trie<String> codeTrie;
	
	public TranslateCharDagVisitor(CharDagVisitor<R> parent, Trie<String> codeTrie, boolean revert) {
		this.parent = parent;
		this.codeTrie = codeTrie;
		this.revert = revert;
	}

	@Override
	public int compareTo(CharDagVisitor<R> o1) {
		TranslateCharDagVisitor<R> o = (TranslateCharDagVisitor<R>)o1;
		int re = Integer.compare(p, o.p);
		if (re==0)
			re = ArrayUtils.compare(codon, o.codon);
		if (re==0)
			re = parent.compareTo(o.parent);
		return re;
	}

	@Override
	public Iterator<R> accept(char c, int pos) {
		codon[p++] = c;
		if (p==3) {
			p=0;
			String aa = codeTrie.get(revert?new CharIterator.ReverseArrayIterator(codon):new CharIterator.ArrayIterator(codon));
			if (aa==null) return parent.accept('X', pos);
			return parent.accept(aa.charAt(0), pos);
		}
		return null;
	}
	
	@Override
	public Iterator<R> accept(char c, SequenceVariant variant, int varToPos) {
		codon[p++] = c;
		if (p==3) {
			p=0;
			String aa = codeTrie.get(revert?new CharIterator.ReverseArrayIterator(codon):new CharIterator.ArrayIterator(codon));
			if (aa==null) return parent.accept('X', variant, varToPos);
			return parent.accept(aa.charAt(0), variant, varToPos);
		}
		return null;
	}
	
//	@Override
//	public void shift(int pos, int shift, SequenceVariant variant) {
//		parent.shift(pos,shift,variant);
//	}

	@Override
	public void addVariant(SequenceVariant variant) {
		parent.addVariant(variant);
	}
	
	@Override
	public CharDagVisitor<R> branch() {
		TranslateCharDagVisitor<R> re = new TranslateCharDagVisitor<>(parent.branch(),codeTrie,revert);
		re.codon = codon.clone();
		re.p = p;
		return re;
	}

	@Override
	public void prune(int pos) {
		parent.prune(pos-p);
	}

	@Override
	public Iterator<SequenceVariant> getVariants() {
		return parent.getVariants();
	}

	

	

}
