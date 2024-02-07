package gedi.util.datastructure.tree.suffixTree.tree;

import java.util.function.Function;
import java.util.function.ToIntFunction;


public class Localization {

	private int node;
	private int infixStart;
	private int infixEnd;

	public Localization(SuffixTree t, int node, int infixStart, int infixEnd) {
		this.node = node;
		this.infixStart = infixStart;
		this.infixEnd = infixEnd;
	}
	
	public int getNode() {
		return node;
	}

	public void setNode(int node) {
		this.node = node;
	}

	public int getInfixStart() {
		return infixStart;
	}

	public void setInfixStart(int infixStart) {
		this.infixStart = infixStart;
	}
	
	public void incrementStart() {
		this.infixStart++;
	}
	
	public void incrementStart(int count) {
		this.infixStart+=count;
	}

	public int getInfixEnd() {
		return infixEnd;
	}

	public void setInfixEnd(int infixEnd) {
		this.infixEnd = infixEnd;
	}
	
	public int getInfixLength() {
		return infixEnd-infixStart;
	}

	
	public boolean isNode() {
		return infixEnd-infixStart==0;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Localization))
			return false;
		Localization o = (Localization) obj;
//		canonize();
//		o.canonize();
		return o.node==node && infixStart==o.infixStart && infixEnd==o.infixEnd;
	}
	
	@Override
	public int hashCode() {
//		canonize();
		return (node|(1<<12))+(infixStart|(1<<6))+infixEnd;
	}
	
	@Override
	public String toString() {
		return "("+node+",\""+infixStart+" "+infixEnd+"\")";
	}



	public static class ToNodeTransformer implements ToIntFunction<Localization> {

		@Override
		public int applyAsInt(Localization arg0) {
			return arg0.getNode();
		}
		
	}

	public void removeInfix() {
		infixStart = 0;
		infixEnd = 0;
	}
	
	
}
