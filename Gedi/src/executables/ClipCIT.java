package executables;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.function.UnaryOperator;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.AlignedReadsDeletion;
import gedi.core.data.reads.AlignedReadsInsertion;
import gedi.core.data.reads.AlignedReadsMismatch;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.StringUtils;
import gedi.util.functions.ExtendedIterator;
import gedi.util.userInteraction.progress.ConsoleProgress;


public class ClipCIT {

	private static int checkIntParam(String[] args, int index) {
		String re = checkParam(args, index);
		if (!StringUtils.isInt(re)) throw new RuntimeException("Must be an integer: "+args[index-1]);
		return Integer.parseInt(args[index]);
	}
	
	private static String checkParam(String[] args, int index)  {
		if (index>=args.length || args[index].startsWith("-")) throw new RuntimeException("Missing argument for "+args[index-1]);
		return args[index];
	}
	
	private static GenomicRegion checkRange(String[] args, int index)  {
		String r = checkParam(args, index);
		if (r.contains("-")) return GenomicRegion.parse(r);
		return new ArrayGenomicRegion(checkIntParam(args, index),Integer.MAX_VALUE);
	}
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws IOException {
		if (args.length<3) {
			usage();
			System.exit(1);
		}
		
		boolean progress = false;
		String out = null;
		GenomicRegion first = null;
		
		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-p"))
				progress = true;
			else if (args[i].equals("-c")) {
				first = checkRange(args, ++i);
			} 
			else {
				out = args[i++];
				args = Arrays.copyOfRange(args, i, args.length);
				i = args.length;
			}
		}
		if (out==null) {
			usage();
			System.exit(1);
		}
		
		if (args.length!=1) {
			usage();
			System.exit(1);
		}
			
		
		if (new File(out).equals(new File(args[0]))) {
			usage();
			System.exit(1);
		}
		
		CenteredDiskIntervalTreeStorage<AlignedReadsData> in =(CenteredDiskIntervalTreeStorage.load(args[0])); 
		CenteredDiskIntervalTreeStorage<AlignedReadsData> outCit = new CenteredDiskIntervalTreeStorage<>(out,in.getType(),in.isCompressed());
		outCit.setMetaData(in.getMetaData());
		
		AlignedReadClipperStartToEnd clipper = new AlignedReadClipperStartToEnd(first.getStart(), first.getEnd());

		ExtendedIterator<? extends ImmutableReferenceGenomicRegion<AlignedReadsData>> it = in.ei();
		if (progress)
			it = it.progress(new ConsoleProgress(System.err), (int)in.size(), a->a.toLocationString());
		it = it.map(a->clipper.apply(a));
		
		outCit.fill(it);
		
	}

	private static void usage() {
		System.out.println("ClipCIT [-p] [-c start[-end]] <output> <file1> \n\n -c is for both first and second read\n -p shows progress");
	}
	
	

	public static class AlignedReadClipperStartToEnd implements UnaryOperator<ImmutableReferenceGenomicRegion<AlignedReadsData>> {
	
		private ArrayGenomicRegion clip;
		
		public AlignedReadClipperStartToEnd(int start, int end) {
			clip = new ArrayGenomicRegion(start,end);
		}
		
		
		
		@Override
		public ImmutableReferenceGenomicRegion<AlignedReadsData> apply(
				ImmutableReferenceGenomicRegion<AlignedReadsData> t) {
			
			ArrayGenomicRegion cclip = clip.intersect(new ArrayGenomicRegion(0,t.getRegion().getTotalLength()));
			ArrayGenomicRegion reg = t.getRegion().intersect(t.map(cclip));
			
			AlignedReadsDataFactory fac = new AlignedReadsDataFactory(t.getData().getNumConditions()).start();
			for (int d=0; d<t.getData().getDistinctSequences(); d++)
				fac.add(t.getData(), d, r->transformVariation(r,cclip));
			
			fac.makeDistinct();
			ImmutableReferenceGenomicRegion<AlignedReadsData> re = new ImmutableReferenceGenomicRegion<>(t.getReference(), reg, fac.create());
			checkCount(t.getData(),re.getData());
			
			return re;
		}
		
		public AlignedReadsVariation transformVariation(AlignedReadsVariation var, ArrayGenomicRegion clip) {
			if (var.isSoftclip()) {
				return null;
			}
			
			
			if (var.getPosition()<clip.getStart() || var.getPosition()>=clip.getEnd())
				return null;
			
			if (var.isMismatch()) {
				if (clip.getStart()>0)
					return new AlignedReadsMismatch(var.getPosition()-clip.getStart(), var.getReferenceSequence(), var.getReadSequence(), var.isFromSecondRead());
				return var;
			}
			else if(var.isDeletion()) {
				if (clip.getStart()>0)
					return new AlignedReadsDeletion(var.getPosition()-clip.getStart(), var.getReferenceSequence(), var.isFromSecondRead());
				return var;
			}
			else if(var.isInsertion()) {
				if (clip.getStart()>0)
					return new AlignedReadsInsertion(var.getPosition()-clip.getStart(), var.getReadSequence(), var.isFromSecondRead());
				return var;
			}
			throw new RuntimeException("Variant not supported!");
		}
	
	
		private void checkCount(AlignedReadsData a, AlignedReadsData b) {
			if (a.getTotalCountOverallInt(ReadCountMode.All)!=b.getTotalCountOverall(ReadCountMode.All))
				throw new RuntimeException();
		}
	
	
	}
}
