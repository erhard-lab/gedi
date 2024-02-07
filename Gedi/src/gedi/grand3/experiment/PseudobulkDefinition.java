package gedi.grand3.experiment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;

import gedi.util.ArrayUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.math.stat.counting.Counter;

public class PseudobulkDefinition {

	private String[] pseudobulkNames;
	private int[][] indicesToPseudobulk;
	private int[] pseudoBulkToSample;
	
	
	public PseudobulkDefinition(String pseudobulkFile, ExperimentalDesign design, Logger log, double minimalPurity) throws IOException {
		HashMap<String, Integer> nameToIndex = EI.seq(0, design.getCount()).map(i->design.getFullName(i)).indexPosition();
		HashMap<String,Integer> newNameToIndex = new HashMap<>();
		IntArrayList[] mapConstruct = new IntArrayList[design.getCount()];
		for (int i=0; i<mapConstruct.length; i++) mapConstruct[i] = new IntArrayList();
		ArrayList<IntArrayList> conditionToSample = new ArrayList<IntArrayList>(); // this stores the sample id for the new condition id; i.e. this must be unique!
		
		HeaderLine h = new HeaderLine();
		for (String[] cols : EI.lines(pseudobulkFile).header(h,'\t',"Cell","Pseudobulk").split('\t').loop()) {
			Integer id=nameToIndex.get(cols[h.get("Cell")]);
			if (id==null) throw new RuntimeException("Cell named "+cols[h.get("Cell")]+" unknown!");
			int mapped = newNameToIndex.computeIfAbsent(cols[h.get("Pseudobulk")],x->newNameToIndex.size());
			if (mapped==conditionToSample.size()) {
				conditionToSample.add(new IntArrayList());
			}
			
			mapConstruct[id].add(mapped);
			conditionToSample.get(mapped).add(design.getIndexToSampleId()[id]);
		}
		pseudobulkNames = new String[newNameToIndex.size()];
		for (String n : newNameToIndex.keySet())
			pseudobulkNames[newNameToIndex.get(n)]=n;
		
		pseudoBulkToSample = new int[newNameToIndex.size()];
		IntArrayList remove = new IntArrayList();
		for (int i=0; i<newNameToIndex.size(); i++) {
			// all cells for a mapped condition have to belong to a unique sample (otherwise estimation is not possible at the moment)
			// if 100% belong to a single sample, fine
			// if it's at least 95%, issue a warning and proceed anyways
			// exclude otherwise!
			Counter<Integer> tab = EI.wrap(conditionToSample.get(i)).tabulate();
			if (tab.elements().size()==1) {} // 100% -> fine
			else if (tab.get(tab.getMaxElement(0),0)>=minimalPurity*tab.total(0)) {
				if (log!=null) log.warning("There are "+String.format("%.0f", (1-tab.get(tab.getMaxElement(0),0)/(double)tab.total(0))*100)+"% pseudobulk cells from other samples for "+pseudobulkNames[i]+". Proceeding anyways (as it is > "+String.format("%.0f", (minimalPurity)*100)+"%)!");
			} else {
				if (log!=null) log.severe("There are "+String.format("%.0f", (1-tab.get(tab.getMaxElement(0),0)/(double)tab.total(0))*100)+"% pseudobulk cells from other samples for "+pseudobulkNames[i]+". Removing (as it is < "+String.format("%.0f", (minimalPurity)*100)+"%)!");
				remove.add(i);
			}
			pseudoBulkToSample[i] = tab.getMaxElement(0);
		}
		
		for (int ri=remove.size()-1; ri>=0; ri--) {
			pseudobulkNames = ArrayUtils.removeIndexFromArray(pseudobulkNames, remove.getInt(ri));
			pseudoBulkToSample = ArrayUtils.removeIndexFromArray(pseudoBulkToSample, remove.getInt(ri));
		}
		
		indicesToPseudobulk = new int[design.getCount()][];
		int[] reindex = new int[newNameToIndex.size()];
		int ind = 0;
		for (int i=0; i<reindex.length; i++) {
			if (!remove.contains(i))
				reindex[i] = ind++;
			else
				reindex[i] = -1;
		}
		if (ind!=pseudobulkNames.length) throw new RuntimeException("Assertion failed!");
		
		for (int i=0; i<indicesToPseudobulk.length; i++) {
			indicesToPseudobulk[i] = mapConstruct[i].iterator()
					.mapIntToInt(inde->reindex[inde])
					.filterInt(ii->ii>=0)
					.toIntArray();
		}
		
	}
	
	public int getNumPseudobulks() {
		return pseudobulkNames.length;
	}
	
	public String getPseudobulkName(int p) {
		return pseudobulkNames[p];
	}
	
	public String[] getPseudobulkNames() {
		return pseudobulkNames;
	}
	
	public int getSampleForPseudobulk(int p) {
		return pseudoBulkToSample[p];
	}
	
	public int[] getSampleForPseudobulks() {
		return pseudoBulkToSample;
	}
	
	public int[] getCellsForPseudobulk(int p) {
		return indicesToPseudobulk[p];
	}
	
	public int[][] getCellsToPseudobulk() {
		return indicesToPseudobulk;
	}
	

}
