package gedi.grand3.targets;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

import gedi.core.reference.ReferenceSequence;
import gedi.util.StringUtils;
import gedi.util.functions.EI;

public class SnpData {

	private final HashMap<String, HashSet<Integer>> data;
	private int size;
	
	public SnpData(File file) throws IOException {
		data = new HashMap<String, HashSet<Integer>>();
		for (String loc : EI.lines(file).splitField('\t', 0).skip(1).loop()) {
			String name = StringUtils.splitField(loc, ':', 0);
			int pos = Integer.parseInt(StringUtils.splitField(loc, ':', 1));
			data.computeIfAbsent(name, x->new HashSet<>()).add(pos);
			size++;
		}
	}
	
	public int size() {
		return size;
	}

	public boolean isSnp(ReferenceSequence refInd, int pos) {
		HashSet<Integer> set = data.get(refInd.getName());
		return set!=null && set.contains(pos);
	}
	
	
	
}
