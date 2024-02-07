package executables;

import java.io.File;
import java.io.IOException;

import javax.script.ScriptException;

import org.broad.igv.bbfile.BBFileReader;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;

import gedi.core.data.numeric.diskrmq.DiskGenomicNumericBuilder;
import gedi.core.reference.Chromosome;
import gedi.core.reference.Strand;
import gedi.util.StringUtils;
import gedi.util.nashorn.JSFunction;
import gedi.util.userInteraction.progress.ConsoleProgress;

public class BW2RMQ {

	public static void main(String[] args) throws IOException, ScriptException {
		if (args.length<2) {
			usage();
			System.exit(1);
		}
		
		if (!new File(args[0]).exists()) {
			System.out.println("File "+args[0]+" does not exist!");
			usage();
			System.exit(1);
		}
		
		BBFileReader[] wigs = {new BBFileReader(args[0])};
		if (args[0].contains("forward") && new File(args[0].replace("forward", "reverse")).exists())
			wigs = new BBFileReader[] {new BBFileReader(args[0]), new BBFileReader(args[0].replace("forward", "reverse"))};
		JSFunction fun = null;
		
		if (args.length==3) 
			fun = new JSFunction(args[1]);
		
		DiskGenomicNumericBuilder rmq = new DiskGenomicNumericBuilder(args[args.length-1],false);
		rmq.setReferenceSorted(true);
		
		ConsoleProgress pr = new ConsoleProgress();
		pr.init();

		int valsadded = 0;
		for (int i=0; i<wigs.length; i++) {
			BBFileReader wig = wigs[i];
			for (String ref : wig.getChromosomeNames()) {
				pr.setDescription(ref+" values: "+valsadded);
				Chromosome rs = Chromosome.obtain(StringUtils.removeFooter(ref,".1")+Strand.values()[i]);
				int id = wig.getChromosomeID(ref);
				BigWigIterator it = wig.getBigWigIterator(ref, wig.getChromosomeBounds(id, id).getStartBase(), ref, wig.getChromosomeBounds(id, id).getEndBase(), false);
				while (it.hasNext()) {
					WigItem item = it.next();
					for (int p=item.getStartBase(); p<item.getEndBase(); p++) {
						float v = item.getWigValue();
						if (v!=0) {
							if (fun==null)
								rmq.addValue(rs, p, v);
							else
								rmq.addValue(rs, p, (Number)fun.apply(v));
							valsadded++;
						}
					}
					
					pr.incrementProgress();
				}
			}
		}
		
		pr.finish();		
		
		rmq.build(false, true);
		
	}

	private static void usage() {
		System.out.println("BW2RMQ <bw> <js-function> <rmq>");
	}
	
}
