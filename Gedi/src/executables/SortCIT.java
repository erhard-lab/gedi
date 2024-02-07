package executables;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.FileUtils;
import gedi.util.functions.EI;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;

/**
 * Most basic by loading it chromosome by chromosome into memory (could be done by repeated intersection queries!)
 * @author erhard
 *
 */
public class SortCIT {

	
	public static void main(String[] args) throws IOException {
		boolean progress = false;
		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else {
				args = Arrays.copyOfRange(args, i, args.length);
				i = args.length;
			}
		}
		if (args.length!=1) {
			usage();
			System.exit(1);
		}
		
		String outPath = FileUtils.insertSuffixBeforeExtension(args[args.length-1],".atac");
		
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> in = new CenteredDiskIntervalTreeStorage<>(args[args.length-1]);
		CenteredDiskIntervalTreeStorage<DefaultAlignedReadsData> out = new CenteredDiskIntervalTreeStorage<>(outPath,in.getType(),in.isCompressed());

		ArrayList<ReferenceSequence> refs = new ArrayList<>(in.getReferenceSequences());
		Collections.sort(refs);
		
		boolean uprogress = progress;
		out.fill(EI.wrap(refs).unfold(ref->{
			ArrayList<ImmutableReferenceGenomicRegion<DefaultAlignedReadsData>> l = in.ei(ref).progress(uprogress?new ConsoleProgress(System.err):new NoProgress(),(int)in.size(ref), d->d.toString()).list();
			Collections.sort(l);
			return EI.wrap(l);
		}));
		
	}


	private static void usage() {
		System.out.println("SortCIT [-p] <cit> ... \n\n -p shows progress");
	}
	
}
