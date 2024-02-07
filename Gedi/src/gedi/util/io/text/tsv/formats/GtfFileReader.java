package gedi.util.io.text.tsv.formats;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.ToIntBiFunction;

import gedi.core.data.annotation.Transcript;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.GeneralUtils;
import gedi.util.StringUtils;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.tsv.GenomicExonsTsvFileReader;
import gedi.util.parsing.IntegerParser;

public class GtfFileReader extends GenomicExonsTsvFileReader<Transcript> {

	
	
	

	public GtfFileReader(String path, String... gtfFeatureType) throws IOException {
		super(path,false,"\t",getRef(),getStart(),getEnd(),getData(),getComp(),init(), Transcript.class);
			
		idgetter = d->d.getGeneId()+","+d.getTranscriptId();
//		readCoordinates= (h,f)->f[2].equals("exon");
		if (gtfFeatureType.length>0) {
			HashSet<String> features = new HashSet<String>(Arrays.asList(gtfFeatureType));
			readCoordinates= (h,f)->getGtfField("transcript_id",f)!=null && features.contains(f[2]);
		} else
			readCoordinates= (h,f)->getGtfField("transcript_id",f)!=null;
		mergeOverlap = true;
		
		lineChecker = check().fieldCount(9)
				.fieldType(3, new IntegerParser())
				.fieldType(4, new IntegerParser())
				.fieldContent(8, f->getGtfField("gene_id",f)==null?"Cannot find gene_id in attributes":null);
//				.fieldContent(8, f->getGtfField("transcript_id",f)==null?"Cannot find transcript_id in attributes":null);
	}
	
	public void setTableOutput(String geneOut, String transOut, String... transColumns) {
		
		LineOrientedFile genes = new LineOrientedFile(geneOut);
		LineOrientedFile trans = new LineOrientedFile(transOut);
		HashSet<String> processedGenes = new HashSet<String>();
		
		init = h->{
				try {
					genes.startWriting();
					genes.writef("Gene ID\tGene Symbol\tBiotype\tSource\n");
					trans.startWriting();
					trans.writef("Transcript ID\tProtein ID\tBiotype\tSource");
					for (String c : transColumns)
						trans.writef("\t%s", c);
					trans.writeLine();
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
		};
		finished = r->{
			try {
				genes.finishWriting();
				trans.finishWriting();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		};
		
		sideEffect = (rgr,h,a)->{
			try {
				if (processedGenes.add(getGtfField("gene_id",a[0]))) 
					genes.writef("%s\t%s\t%s\t%s\n",getGtfField("gene_id",a[0]),getGtfField("gene_name",a[0],""),getGtfField("gene_biotype",a[0],""),getGtfField("gene_source",a[0],""));
				String prot = "";
				String src = getGtfField("transcript_source",a[0],"");
				for (int i=0; i<a.length; i++) {
					if (getGtfField("protein_id",a[i])!=null)
						prot = getGtfField("protein_id",a[i]);
					if (getGtfField("ccds_id",a[i])!=null)
						src = "CCDS"; // technically speaking, this is wrong
				}
				trans.writef("%s\t%s\t%s\t%s",getGtfField("transcript_id",a[0]),prot,getGtfField("transcript_biotype",a[0],""),src);
				for (String c : transColumns)
					trans.writef("\t%s", getGtfField(c,a[0]));
				trans.writeLine();
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		};
	}
	

	private static Consumer<HeaderLine> init() {
		return null;
	}


	private static BiPredicate<String[],String[]> getComp() {
		return (a,b)->GeneralUtils.isEqual(getGtfField("gene_id",a),getGtfField("gene_id",b))&&GeneralUtils.isEqual(getGtfField("transcript_id",a),getGtfField("transcript_id",b));
	}


	private static ToIntBiFunction<HeaderLine, String[]> getEnd() {
		return (h,f)-> {
			return Integer.parseInt(f[4]);
		};
	}


	private static ToIntBiFunction<HeaderLine, String[]> getStart() {
		return (h,f)-> {
			return Integer.parseInt(f[3])-1;
		};
	}


	private static BiFunction<HeaderLine, String[], ReferenceSequence> getRef() {
		return (h,f)-> {
			return Chromosome.obtain(f[0],f[6]);
		};
	}
	
	private static final String getGtfField(String name, String[] f) {
		return getGtfField(name, f[8],null);
	}
	
	public static final String getGtfField(String name, String f8) {
		return getGtfField(name, f8, null);
	}
	
	private static final String getGtfField(String name, String[] f, String def) {
		return getGtfField(name, f[8],def);
	}
	
	public static final String getGtfField(String name, String f8, String def) {
		String[] f = StringUtils.split(f8, ';');
		for (int i=0; i<f.length; i++) {
			if (f[i].contains(name))
				return StringUtils.trim(f[i].substring(f[i].indexOf(name)+name.length()), ' ','"');
		}
		return def;
	}


	private static BiFunction<HeaderLine, String[][], Transcript> getData() {
		return (h,tr)->{
			
			boolean minus = tr[0][6].equals("-");
			int codingStart = -1;
			int codingEnd = -1;
			for (int i=0; i<tr.length; i++) {
				if (tr[i][2].equals("start_codon")) {
					if (!minus) {
						if (codingStart==-1) codingStart = Integer.parseInt(tr[i][3])-1;
						else codingStart = Math.min(codingStart,Integer.parseInt(tr[i][3])-1);
					}
					else {
						if (codingEnd==-1) codingEnd = Integer.parseInt(tr[i][4]);
						else codingEnd = Math.max(codingEnd,Integer.parseInt(tr[i][4]));
					}
				} else if (tr[i][2].equals("stop_codon")) {
					if (!minus) {
						if (codingEnd==-1) codingEnd = Integer.parseInt(tr[i][4]);
						else codingEnd = Math.max(codingEnd,Integer.parseInt(tr[i][4]));
					}
					else {
						if (codingStart==-1) codingStart = Integer.parseInt(tr[i][3])-1;
						else codingStart = Math.min(codingStart,Integer.parseInt(tr[i][3])-1);
					}
				} else if (tr[i][2].equals("CDS")) {
					if (codingStart==-1) codingStart = Integer.parseInt(tr[i][3])-1;
					else codingStart = Math.min(codingStart,Integer.parseInt(tr[i][3])-1);
					if (codingEnd==-1) codingEnd = Integer.parseInt(tr[i][4]);
					else codingEnd = Math.max(codingEnd,Integer.parseInt(tr[i][4]));
				}
			}
			Transcript re = new Transcript(
					getGtfField("gene_id", tr[0]),
					getGtfField("transcript_id", tr[0]),
					codingStart,codingEnd
					);
			return re;
		};
	}

	public static void main(String[] args) throws IOException {
		for (ImmutableReferenceGenomicRegion<Transcript> r :new GtfFileReader(args[0]).readIntoMemoryThrowOnNonUnique().getReferenceGenomicRegions() )
			System.out.println(r);
		
		
	}
}
