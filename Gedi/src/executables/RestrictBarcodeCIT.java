package executables;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.util.ParseUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.mutable.MutableInteger;
import gedi.util.sequence.DnaSequence;
import gedi.util.userInteraction.progress.ConsoleProgress;

import java.io.File;
import java.io.IOException;

public class RestrictBarcodeCIT {

	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		
		boolean progress = false;
		
		int i;
		for (i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else
				break;
		}
		
		if (i+3!=args.length) {
			usage();
			System.exit(1);
		}
		
		CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData> in = new CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData>(args[i++]);
		
		int numCond = in.getRandomRecord().getNumConditions();
		if (numCond!=1) throw new RuntimeException("Restricting barcodes is not implemented for mult-condition CIT files!");
		
		
		Trie<Boolean> barcodes = new Trie<Boolean>();
		for (String bc : (new File(args[i]).exists()? EI.lines(args[i]) : EI.split(args[i], ',')).loop())
			barcodes.put(bc, true);
		i++;
		
		System.err.println("Filtering for "+barcodes.size()+" barcodes!");
				
		CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<BarcodedAlignedReadsData>(args[i++],BarcodedAlignedReadsData.class);
		
		
		out.fill(in.ei()
				.iff(progress, e->e.progress(new ConsoleProgress(System.err),(int)in.size(),rgr->rgr.toLocationStringRemovedIntrons()))
				.map(rgr-> {
					return rgr.toMutable().transformData(ard->restrict(ard,barcodes));
				})
				.filter(rgr->rgr.getData()!=null),
				new ConsoleProgress(System.err)
				);
		
		
		
		StringBuilder sb = new StringBuilder();
		sb.append("{");
		for (String s : in.getMetaData().getProperties()) {
			if (s.equals("name")) {
				if (sb.length()>1) sb.append(",");
				
				sb.append("\"").append(s).append("\":");
				DynamicObject e = in.getMetaData().getEntry(s);
				sb.append(in.getMetaData().getEntry(s).toJson());
			}
		}
		sb.append("}");
		out.setMetaData(DynamicObject.parseJson(sb.toString()));
		
		
	}

	private static BarcodedAlignedReadsData restrict(BarcodedAlignedReadsData ard, Trie<Boolean> barcodes) {
		for (int d=0; d<ard.getDistinctSequences(); d++) {
			for (DnaSequence bcx : ard.getBarcodes(d, 0)) {
				if (barcodes.iteratePrefixMatches(bcx).count()==0) {
					// ok, there is a sequence that I have to remove, so do the factory thing!
					AlignedReadsDataFactory fac = new AlignedReadsDataFactory(1);
					fac.start();
					for (d=0; d<ard.getDistinctSequences(); d++) {
						DnaSequence[] bcs = EI.wrap(ard.getBarcodes(d, 0)).filter(bcc->barcodes.iteratePrefixMatches(bcc).count()>0).toArray(DnaSequence.class);
						if (bcs.length>0) {
							fac.add(ard, d);
							fac.setCount(fac.getDistinctSequences()-1,0, bcs.length, bcs);
						}
					}
					
					return fac.getDistinctSequences()==0?null:fac.createBarcode();
				}
			}
		}
		
		
		// if we are here, then all barcodes in ard are valid, so nothing to do!
		return ard;
	}

	private static void usage() {
		System.out.println("RestrictBarcodeCIT [-p] <input> <barcodes> <output> \n\n -p shows progress\n barcodes either is a comma-separated list of barcodes, or a file that line by line contains barcodes to keep. Barcodes must have the same length!");
	}
	
}
