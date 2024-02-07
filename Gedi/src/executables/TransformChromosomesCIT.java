package executables;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.reference.Chromosome;
import gedi.core.reference.Strand;
import gedi.util.StringUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class TransformChromosomesCIT {

	public static void main(String[] args) throws IOException {
		if (args.length<1) {
			usage();
			System.exit(1);
		}
		
		if (!new File(args[0]).exists()) {
			System.out.println("File "+args[0]+" does not exist!");
			usage();
			System.exit(1);
		}
		
		if (args.length<2) {
			System.out.println("No mapping given!");
			usage();
			System.exit(1);
		}
		boolean reverse = args[1].equals("-r");
		int start = reverse?2:1;
		
		HashMap<String,String> mapping = new HashMap<String, String>();
		for (int i=start; i<args.length; i++) {
			String[] p = StringUtils.split(args[i], "->");
			if (p.length!=2) {
				System.out.println("Wrong mapping format: "+args[i]);
				usage();
				System.exit(1);
			}
			mapping.put(p[0], p[1]);
		}
		
		
		new CenteredDiskIntervalTreeStorage<Object>(args[0]).mapChromosomes(ref->{
			Strand strand = reverse?ref.getStrand().toOpposite():ref.getStrand();
			return Chromosome.obtain(mapping.containsKey(ref.getName())?mapping.get(ref.getName()):ref.getName(),strand);
		});
		
		
	}

	private static void usage() {
		System.out.println("TransformChromosomesCIT <file> [-r] [oldname->newname ...]");
	}
	
}
