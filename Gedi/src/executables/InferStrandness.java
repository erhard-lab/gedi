package executables;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;

import gedi.app.Gedi;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.GenomicRegionStoragePreload;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.core.workspace.loader.WorkspaceItemLoader;
import gedi.core.workspace.loader.WorkspaceItemLoaderExtensionPoint;
import gedi.util.StringUtils;


public class InferStrandness {

	private static String checkParam(String[] args, int index) {
		if (index>=args.length || args[index].startsWith("-")) throw new RuntimeException("Missing argument for "+args[index-1]);
		return args[index];
	}
	private static int checkMultiParam(String[] args, int index, int max, ArrayList<String> re) {
		while (index<max && !args[index].startsWith("-")) 
			re.add(args[index++]);
		return index-1;
	}

	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		
		Gedi.startup(false);
		boolean full = false;
		Genomic g = null;
		
		int i;
		for (i=0; i<args.length; i++) {
			if (args[i].equals("-h")) {
				usage();
				return;
			}
			else if (args[i].equals("-g")) {
				ArrayList<String> gnames = new ArrayList<>();
				i = checkMultiParam(args, ++i, args.length-1,gnames);
				g = Genomic.get(gnames);
			}
			else if (args[i].equals("-full"))
				full = true;
			else
				break;
		}
		
		if (i+1!=args.length) {
			usage();
			System.exit(1);
		}
		
		
		Path p = Paths.get(args[i++]);
		WorkspaceItemLoader<GenomicRegionStorage<? extends AlignedReadsData>,GenomicRegionStoragePreload<AlignedReadsData>> loader = WorkspaceItemLoaderExtensionPoint.getInstance().get(p);
		if (loader==null)
			throw new RuntimeException("No loader available for "+p);
		GenomicRegionStorage<? extends AlignedReadsData> in = (GenomicRegionStorage<? extends AlignedReadsData>) loader.load(p);
		
		MemoryIntervalTreeStorage<String> genes = g.getGenes();
		
		int s = 0;
		int as = 0;
		for ( ImmutableReferenceGenomicRegion<? extends AlignedReadsData> r : in.ei().iff(!full,ei->ei.head(500000)).loop()) {
				long sense = genes.ei(r).count();
				long antisense = genes.ei(r.toOppositeStrand()).count();
				if (sense>0 && antisense==0) s++;
				if (sense==0 && antisense>0) as++;
		}
		
		if (s/2>as) 
			System.out.println(Strandness.Sense);
		else if (as/2>s)
			System.out.println(Strandness.Antisense);
		else 
			System.out.println(Strandness.Unspecific);
		
	}

	private static void usage() {
		System.out.println("ReadCount [-g <genomes>] [-p] [-nometa] [-m <mode>]<input>\n\n -mito don't skip mitochondrial reads\n -g count genome specific reads as well\n -p shows progress\n -m <mode> set read count mode  (One of: "+StringUtils.concat(',', ReadCountMode.values())+")\n\nOutputs a table and writes the total counts into the metadata file");
	}
	
}
