package executables;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.util.StringUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

public class CompressCIT {

	public static void main(String[] args) throws IOException {
		if (args.length!=2) {
			usage();
			System.exit(1);
		}
		
		if (!new File(args[0]).exists()) {
			System.out.println("File "+args[0]+" does not exist!");
			usage();
			System.exit(1);
		}
		
		
		CenteredDiskIntervalTreeStorage in = new CenteredDiskIntervalTreeStorage(args[0]);
		CenteredDiskIntervalTreeStorage out = new CenteredDiskIntervalTreeStorage(args[1],in.getType(),true);
		out.fill(in);
		out.setMetaData(in.getMetaData());
		
	}

	private static void usage() {
		System.out.println("CompressCIT <in> <out>");
	}
	
}
