package gedi.core.data.reads.subreads;

import java.io.File;
import java.util.ArrayList;
import java.util.logging.Logger;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.FileUtils;
import gedi.util.functions.BiIntConsumer;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.HeaderLine;



/**
 * This provides a general way to convert all reads in a certain genomic region to subread. For single end and paired end reads, this
 * is unnecessary (this is why {@link SingleEndToSubreadsConverter} and {@link PairedEndToSubreadsConverter} inherits from {@link ReadByReadToSubreadsConverter}
 * @author erhard
 *
 */
public interface ToSubreadsConverter<A extends AlignedReadsData> {

	ToSubreadsConverter<A> setDebug(boolean debug);
	String[] getSemantic();
	boolean isReadByRead();
	
	default void logUsedTotal(Logger logger, int used, int total) {
		logger.info(String.format("Reads or UMIs used/encountered: %d/%d (%.2f%%)\n", used, total,used*100.0/total));
		if (total==0) logger.severe("No reads or UMIs encountered at all!");
	}
	
	default ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> convert(
			ImmutableReferenceGenomicRegion<? extends A> read, boolean sense, MismatchReporter reporter,
			BiIntConsumer usedTotal) {
		throw new RuntimeException("can't do!");
	}

	default ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> convert(
			String id,
			ReferenceSequence reference, 
			ArrayList<ImmutableReferenceGenomicRegion<A>> cluster,
			MismatchReporter reporter,
			BiIntConsumer usedTotal) {
		return EI.wrap(cluster).
				map(r->convert(r,r.getReference().getStrand().equals(reference.getStrand()),reporter,usedTotal));
	}
	
	public static <A extends AlignedReadsData> ToSubreadsConverter<A> infer(GenomicRegionStorage<A> reads, boolean debug) {
		if (reads.getRandomRecord() instanceof SubreadsAlignedReadsData) {
			String[] subreads = reads.getMetaData().getEntry("subreads").asArray(String.class, d->d.asString());
			return (ToSubreadsConverter<A>) new NoopToSubreadsConverter(subreads);
		}
		if (reads.getRandomRecord() instanceof BarcodedAlignedReadsData) {
			String bcFile = FileUtils.getExtensionSibling(reads.getPath(), "barcodes.tsv");
			HeaderLine h = new HeaderLine();
			String[] conditions = new File(bcFile).exists()?EI.lines2(bcFile).header(h).split('\t').map(a->a[h.get("Library")]).toArray(String.class):null;
			String[] cells = new File(bcFile).exists()?EI.lines2(bcFile).header(h).split('\t').map(a->a[h.get("Barcode")]).toArray(String.class):null;
			
			if (cells==null) {
				throw new RuntimeException("Site specific UMI handling not implemented yet!");
//				SiteSpecificSingleEndUmiToSubreadsConverter re = new SiteSpecificSingleEndUmiToSubreadsConverter();
//				re.setDebug(debug);
//				return (ToSubreadsConverter<A>) re;
			}
			
			GeneSpecificUmiSenseToSubreadsConverter re = new GeneSpecificUmiSenseToSubreadsConverter(reads.getMetaDataConditions(),conditions,cells);
			re.setDebug(debug);
			return (ToSubreadsConverter<A>) re;
	
		}
		else if (reads.getRandomRecord().hasGeometry()) {
			PairedEndToSubreadsConverter re = new PairedEndToSubreadsConverter(true);
			re.setDebug(debug);
			return (ToSubreadsConverter<A>) re;
		}
		else {
			SingleEndToSubreadsConverter re = new SingleEndToSubreadsConverter();
			re.setDebug(debug);
			return (ToSubreadsConverter<A>) re;
		}
	}
	
	
	
}