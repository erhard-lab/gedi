package gedi.util.io.text.tsv.formats;

import java.io.IOException;
import java.net.URL;
import java.util.HashSet;
import java.util.Iterator;
import java.util.function.Function;

import gedi.core.data.annotation.Gff3Element;
import gedi.core.data.annotation.NameAnnotation;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.functions.EI;

public class MirbaseGffReader extends GffFileReader<NameAnnotation> {

	public MirbaseGffReader(String path, boolean pre, boolean mature) {
		super(path, getMapper(pre,mature), NameAnnotation.class);
	}

	private static Function<Gff3Element, NameAnnotation> getMapper(boolean pre, boolean mature) {
		HashSet<String> types = new HashSet<String>();
		if (pre) types.add("miRNA_primary_transcript");
		if (mature) types.add("miRNA");
		return (gff)->types.contains(gff.getFeature())?new NameAnnotation((String)gff.getAttribute("Name")):null;
	}
	
	public static GenomicRegionStorage<NameAnnotation> fromWebMature(String org) throws IOException {
		return new MirbaseGffReader("MIRBASE", false, true) {
			@Override
			protected Iterator<String> createIterator() throws IOException {
				return EI.lines(new URL("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/"+org+".gff3").openStream()).filter(s->!s.startsWith("#"));
			}
		}.readIntoMemoryThrowOnNonUnique();
	}

	
	public static GenomicRegionStorage<NameAnnotation> fromWebPrecursor(String org) throws IOException {
		return new MirbaseGffReader("MIRBASE", true, false) {
			@Override
			protected Iterator<String> createIterator() throws IOException {
				return EI.lines(new URL("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/"+org+".gff3").openStream()).filter(s->!s.startsWith("#"));
			}
		}.readIntoMemoryThrowOnNonUnique();
	}
	
	public static GenomicRegionStorage<NameAnnotation> fromFileMature(String... files) throws IOException {
		GenomicRegionStorage<NameAnnotation> re = new MemoryIntervalTreeStorage<NameAnnotation>(NameAnnotation.class);
		
		for (String file : files)
			re.fill(new MirbaseGffReader(file, false, true).readIntoMemoryThrowOnNonUnique());
		
		return re;
	}

	
	public static GenomicRegionStorage<NameAnnotation> fromFilePrecursor(String file) throws IOException {
		return new MirbaseGffReader(file, true, false).readIntoMemoryThrowOnNonUnique();
	}

}
