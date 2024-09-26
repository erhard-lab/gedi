package gedi.grand3.experiment;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gedi.core.data.reads.AlignedReadsData;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;



/**
 * Grand3 is called with a single cit file, call it the experimental design or the experiment. This can contain severak libraries (in cit they are called conditions).
 * 
 * Each library might have several barcodes (e.g. 10x), 
 * barcodes could be partitioned into samples (multi-seq).
 * 
 * This is a tree with four levels: 
 * Experiment
 *   library 1
 *     sample 1.1
 *       barcode 1.1.1
 *       ...
 *       barcode 1.1.k
 *     sample 1.2L
 *     ...
 *     sample 1.l
 *   library 2
 *   ...
 *   
 *   
 * Barcode and sample names but not library names can be empty strings. 
 *   
 * A bulk experiment (or SMART-seq2 based scSLAM) usually has several libraries, each with a single sample and a single barcode (at least logically)
 * A simple 10x experiment has one library, one sample and thousands of barcodes
 * A multi-seq 10x experiment has one library, several samples (multi-seq barcodes) with dozens to hundreds of barcodes each.
 * Each 10x experiment also might involve multiple libraries.
 * 
 * Global parameter are inferred on the level of samples. Target specific NTRs on the level of barcodes.
 * 
 * Each sample has one or several metabolic labels (i.e. is this a no4sU or 4sU library, maybe also 6sG or both).
 * 
 * There are two ways how conditions (as defined by {@link AlignedReadsData}) are numbered:
 * Before conversion as in the original cit file, and after conversion (e.g. barcode expansion).
 * The former is called library id, the latter is called index.
 * 
 * This class also defines a numbering of the samples (across libraries). This is called sample id.
 * 
 * @author flo
 *
 */
public class ExperimentalDesign {

	private MetabolicLabelType[] typeOrder;

	private HashMap<MetabolicLabelType,IntArrayList> labelNameToLibraryId = new HashMap<>();
	private HashMap<MetabolicLabelType,IntArrayList> labelNameToSamples = new HashMap<>();
	private HashMap<MetabolicLabelType,IntArrayList> labelNameToSamplesNot = new HashMap<>();
	
	private String[] libraries;
	private String[] samples;
	private String[] barcodes;
	private MetabolicLabel[][] labels;
	
	private int[] indexToIndex; // a dummy mapper for using it just like indexToSampleId
	private int[] indexToSampleId;
	private int[] indexToLibraryId;
	private int[][] sampleIdToIndices;
	private int[][] libraryIdToIndices;
	
	private MetabolicLabel[][] labelsOfSamples; // as many as there are samples!
	
	
	public ExperimentalDesign(String[] libraries, String[] samples, String[] barcodes, MetabolicLabel[][] labels) {
		
		int n = libraries.length;
		
		if (n!=samples.length || n!=labels.length || n!=barcodes.length)
			throw new RuntimeException("Fatal error in experimental design!");
		EI.wrap(labels).mapToInt(a->a.length).unique(false).getUnique("Error in label definition", "Error in label definition!");
		
		
		this.libraries = libraries;
		this.samples = samples;
		this.barcodes = barcodes;
		this.labels = labels;
		
		this.typeOrder = new MetabolicLabelType[labels[0].length];
		for (int i=0; i<n; i++) 
			for (int l=0; l<labels[i].length; l++) 
				if (labels[i][l]!=null)
					typeOrder[l]=labels[i][l].getType();
		
		for (int l=0; l<typeOrder.length; l++)
			if (typeOrder[l]==null)
				throw new RuntimeException("Metabolic label type undefined; check experimental design!");
		
		this.indexToIndex = new int[n];
		this.indexToSampleId = new int[n];
		this.indexToLibraryId = new int[n];
		
		ArrayList<MetabolicLabel[]> labelsOfSamples = new ArrayList<>();
		
		HashMap<String,Integer> libraryId = new HashMap<String, Integer>();
		HashMap<String,Integer> sampleId = new HashMap<String, Integer>();
		
		boolean noSamples = EI.wrap(samples).filter(s->s.length()==0).count()==samples.length;
		
		boolean firstInSample = false;
		for (int i=0; i<n; i++) {
			firstInSample = false;
			
			indexToIndex[i]=i;
			String sname = libraries[i]+"."+samples[i];
			
			Integer sid = sampleId.get(sname);
			if (sid==null) {
				sampleId.put(sname, sid = sampleId.size());
				firstInSample = true;
				labelsOfSamples.add(labels[i]);
			}
			
			this.indexToSampleId[i] = sid;
			this.indexToLibraryId[i] = libraryId.computeIfAbsent(libraries[i], x->libraryId.size());
			
			// as it is unclear where the labels are defined, just check this here (i.e. in every line)
			for (int l=0; l<labels[i].length; l++) {
				if (labels[i][l]!=null){
					labelNameToLibraryId.computeIfAbsent(labels[i][l].getType(), x->new IntArrayList()).add(libraryId.get(libraries[i]));
					if (firstInSample) 
						labelNameToSamples.computeIfAbsent(labels[i][l].getType(), x->new IntArrayList()).add(sampleId.get(sname));
				} 
				else if (firstInSample) 
					labelNameToSamplesNot.computeIfAbsent(typeOrder[l], x->new IntArrayList()).add(sampleId.get(sname));
				
				// check that labels are homogeneous among samples
				if (!noSamples && labels[i][l]!=labelsOfSamples.get(sampleId.get(sname))[l] && !labels[i][l].equals(labelsOfSamples.get(sampleId.get(sname))[l]))
					throw new AssertionError("Labels are not homogenous in samples!");
			}
		}
		
		for (MetabolicLabelType lab : typeOrder) // if there is no no4sU sample!
			labelNameToSamplesNot.computeIfAbsent(lab, x->new IntArrayList());
		
		for (IntArrayList d : labelNameToLibraryId.values())
			d.unique();
		sampleIdToIndices = index(sampleId.size(), indexToSampleId);
		libraryIdToIndices = index(libraryId.size(), indexToLibraryId);
		
		this.labelsOfSamples = labelsOfSamples.toArray(new MetabolicLabel[0][]);
				
	}

	private int[][] index(int n, int[] a) {
		IntArrayList[] re = new IntArrayList[n];
		for (int i=0; i<re.length; i++)
			re[i] = new IntArrayList(a.length);

		for (int i=0; i<a.length; i++) 
			re[a[i]].add(i);
		
		int[][] re2 = new int[n][];
		for (int i=0; i<re.length; i++)
			re2[i] = re[i].toIntArray();
		
		return re2;
	}

	public MetabolicLabelType[] getTypes() {
		return typeOrder;
	}
	
	public String getFullName(int index) {
		StringBuilder sb = new StringBuilder();
		sb.append(libraries[index]);
		if (samples[index].length()>0) sb.append(".").append(samples[index]);
		if (barcodes[index].length()>0) sb.append(".").append(barcodes[index]);
		return sb.toString();
	}
	
	public String[] getBarcodes() {
		return barcodes;
	}
	
	public String getSampleName(int index) {
		StringBuilder sb = new StringBuilder();
		sb.append(libraries[index]);
		if (samples[index].length()>0) sb.append(".").append(samples[index]);
		return sb.toString();
	}

	public String getSampleNameForSampleIndex(int sample) {
		return getSampleName(sampleIdToIndices[sample][0]);
	}
	
	public int[] getIndexToSampleId() {
		return indexToSampleId;
	}
	public int[] getIndexToIndex() {
		return indexToIndex;
	}
	
	public int[] getIndicesForSampleId(int sample) {
		return sampleIdToIndices[sample];
	}
	
	public int getNumSamples() {
		return sampleIdToIndices.length;
	}
	/**
	 * How many barcodes
	 * @return
	 */
	public int getCount() {
		return libraries.length;
	}


	public int getNumLibraries() {
		return libraryIdToIndices.length;
	}

	
	public MetabolicLabel getLabelForSample(int sampleId, MetabolicLabelType type) {
		int ti = ArrayUtils.linearSearch(typeOrder, type);
		if (ti==-1) throw new RuntimeException("Fatal exception in experimental design: "+type);
		return labelsOfSamples[sampleId][ti];
	}
	public MetabolicLabel getLabelForIndex(int index, MetabolicLabelType type) {
		int ti = ArrayUtils.linearSearch(typeOrder, type);
		if (ti==-1) throw new RuntimeException("Fatal exception in experimental design: "+type);
		return labels[index][ti];
	}


	
	/**
	 * Gets all library ids that have the given metabolic label
	 * @param type
	 * @return
	 */
	public int[] getLibraryIdsWithSamplesHaving(MetabolicLabelType type) {
		return labelNameToLibraryId.get(type).toIntArray();
	}
	
	/**
	 * sample ids!
	 * @param type
	 * @return
	 */
	public int[] getSamplesHaving(MetabolicLabelType type) {
		return labelNameToSamples.get(type).toIntArray();
	}
	/**
	 * sample ids!
	 * @param type
	 * @return
	 */
	public int[] getSamplesNotHaving(MetabolicLabelType type) {
		return labelNameToSamplesNot.get(type).toIntArray();
	}
	
	
	public void writeTable(File f) throws IOException {
		LineWriter out = new LineOrientedFile(f.getPath()).write();
		out.write("Library\tSample\tBarcode");
		for (MetabolicLabelType l : getTypes())
			out.write("\t"+l.toString()+" concentration\t"+l.toString()+" duration\t"+l.toString()+" chase");
		out.writeLine();
		
		int n = libraries.length;
		for (int i=0; i<n; i++) {
			out.writef("%s\t%s\t%s", libraries[i],samples[i],barcodes[i]);
			for (int l=0; l<labels[i].length; l++) {
				if (labels[i][l]==null) out.writef("\tNA\tNA\tNA");
				else out.writef("\t%.0f\t%.2f\t%.2f",labels[i][l].getConcentration(),labels[i][l].getDuration(),labels[i][l].getChase());
			}
			out.writeLine();
		}
		out.close();
	}
	
	public static ExperimentalDesign fromTable(File f) throws IOException {
		ArrayList<String> libraries = new ArrayList<>();
		ArrayList<String> samples = new ArrayList<>();
		ArrayList<String> barcodes = new ArrayList<>();
		ArrayList<MetabolicLabel[]> labels = new ArrayList<>();
		
		ExtendedIterator<String[]> it = EI.lines(f).split('\t');
		String[] header = it.next();
		MetabolicLabelType[] types = EI.wrap(header).filter(s->s.endsWith(" concentration")).map(s->MetabolicLabelType.fromString(s.substring(0,s.indexOf(' ')))).toArray(MetabolicLabelType.class);
		
		if (header[0].equals("Sample")) {
			for (String[] a : it.loop()) {
				libraries.add(a[0]);
				samples.add("");
				barcodes.add(a[1]);
				MetabolicLabel[] label = new MetabolicLabel[types.length];
				for (int i=0; i<types.length; i++)
					if (!a[i*3+2].equals("NA"))
						label[i] = new MetabolicLabel(types[i],Double.parseDouble(a[i*3+2]), Double.parseDouble(a[i*3+3]), Double.parseDouble(a[i*3+4]));
				labels.add(label);
			}
		}
		else {
			for (String[] a : it.loop()) {
				libraries.add(a[0]);
				samples.add(a[1]);
				barcodes.add(a[2]);
				MetabolicLabel[] label = new MetabolicLabel[types.length];
				for (int i=0; i<types.length; i++)
					if (!a[i*3+3].equals("NA"))
						label[i] = new MetabolicLabel(types[i],Double.parseDouble(a[i*3+3]), Double.parseDouble(a[i*3+4]), Double.parseDouble(a[i*3+5]));
				labels.add(label);
			}
		}
		return new ExperimentalDesign(libraries.toArray(new String[0]),samples.toArray(new String[0]),barcodes.toArray(new String[0]),labels.toArray(new MetabolicLabel[0][0]));
	}
	
	
	public static final Pattern regex_no4sU = Pattern.compile("no4sU|no4TU|nos4U|\\.no\\.",Pattern.CASE_INSENSITIVE);
	public static final Pattern[] regex_labels = {
			Pattern.compile("\\.4TU\\.|\\.4sU\\.|\\.s4U\\.",Pattern.CASE_INSENSITIVE),
			Pattern.compile("\\.6TG\\.|6sG\\.|\\.s6G\\.",Pattern.CASE_INSENSITIVE)
	};
	public static final MetabolicLabelType[] types = {MetabolicLabelType._4sU,MetabolicLabelType._6sG};

	public static final Pattern regex_duration = Pattern.compile("([0-9]+)(h|min)(\\.|$)");
	public static final Pattern regex_conc = Pattern.compile("(\\d+)uM(\\.|$)");


	public static ExperimentalDesign infer(Logger log, String[] libraryNames, File barcodeFile) throws IOException {
		
		HashMap<String,String[][]> barcodeMap = new HashMap<String, String[][]>();
		if (barcodeFile!=null) {
			int libraryIndex = 0;
			HeaderLine header = new HeaderLine();
			for (String[][] block : EI.lines(barcodeFile).header(header).split('\t').multiplexUnsorted((a,b)->a[0].equals(b[0]), String[].class).loop()) {
				if (!libraryNames[libraryIndex].equals(block[0][header.get("Library")]))
					log.warning("Library name in barcode file and reads file do not match: "+block[0][header.get("Library")]+" "+libraryNames[libraryIndex]);
				String[] bcs = EI.wrap(block).map(a->a[header.get("Barcode")]).toArray(String.class);
				String[] samples = header.hasField("Sample")?EI.wrap(block).map(a->a[header.get("Sample")]).toArray(String.class):EI.repeat(bcs.length, "").toArray(String.class);
				
				barcodeMap.put(libraryNames[libraryIndex++], new String[][] {samples,bcs});
			}
		}
		else {
			for (String sample:libraryNames)
				barcodeMap.put(sample,new String[][] {new String[] {""},new String[] {""}});
		}
		
		
		int nlabels = EI.wrap(regex_labels).filter(regex->EI.wrap(libraryNames).str().filter(s->!regex_no4sU.matcher(s).find() && regex.matcher(s).find()).count()>0).countInt();
		boolean hasno4sU = false;
		if (nlabels==0)
			for (int cond=0; cond<libraryNames.length; cond++) {
				String[][] two = barcodeMap.get(libraryNames[cond]);
				for (int i=0; i<two[0].length; i++) {
					String c = libraryNames[cond]+"."+two[0][i];
					if (regex_no4sU.matcher(c).find())
						hasno4sU=true;
				}
			}
		
		int total = EI.wrap(barcodeMap.values()).mapToInt(a->a[0].length).sum();
		
//		if (hasno4sU && nlabels>0)
//			throw new RuntimeException("Cannot infer experimental design, there are sample named no4sU AND 4sU/6sG!");
		
		String[] libraries = new String[total];
		String[] samples = new String[total];
		String[] barcodes = new String[total];
		MetabolicLabel[][] labels = new MetabolicLabel[total][];
		int index = 0;
		for (int cond=0; cond<libraryNames.length; cond++) {
			String[][] two = barcodeMap.get(libraryNames[cond]);
			for (int i=0; i<two[0].length; i++) {
				String c = libraryNames[cond]+"."+two[0][i];
				MetabolicLabel[] l = new MetabolicLabel[hasno4sU?1:Math.max(nlabels, 1)];
				
				
				if (hasno4sU) {
					if (!regex_no4sU.matcher(c).find()) {
						l[0]=new MetabolicLabel(types[0],inferNumber(regex_conc,c), inferNumber(regex_duration,c), 0);
						if (l[0].getConcentration()==0 || l[0].getDuration()==0) 
							l[0]=null;
					}
				}
				else if (nlabels==0) {
					l[0]=new MetabolicLabel(types[0],inferNumber(regex_conc,c), inferNumber(regex_duration,c), 0);
					if (l[0].getConcentration()==0 || l[0].getDuration()==0) 
						l[0]=null;
				}
				else {
					for (int r=0; r<regex_labels.length; r++) {
						Matcher m = regex_labels[r].matcher(c);
						if (m.find()) {
								l[r]=new MetabolicLabel(types[r],inferNumber(regex_conc,c), inferNumber(regex_duration,c), 0);
								if (l[r].getConcentration()==0 || l[r].getDuration()==0) 
									l[r]=null;
						} 
					}
				}
			
			
				libraries[index] = libraryNames[cond];
				samples[index] = two[0][i];
				barcodes[index] = two[1][i];
				labels[index++] = l;
			}
				
			
		}
		
		if (EI.wrap(labels).filter(a->a.length>1|a[0]!=null).count()==0) {
			// only no4sU, this would crash grand3, but might be important for just counting reads the grand3 way
			for (int i=0; i<labels.length; i++) 
				labels[i][0] = new MetabolicLabel(types[0],1, 1, 0);
		}
		
		
//		if (
//				EI.wrap(samples).mapToInt(s->s.length()==0?0:1).sum()==0 &&
//						EI.wrap(barcodes).mapToInt(s->s.length()==0?0:1).sum()==0 &&
//						EI.wrap(libraries).mapToInt(s->StringUtils.countChar(s, '.')>=2?0:1).sum()==0) {
//			for (int i=0; i<libraries.length; i++) {
//				int s = libraries[i].indexOf('.');
//				int e = libraries[i].lastIndexOf('.');
//				barcodes[i]=libraries[i].substring(e+1);
//				samples[i]=libraries[i].substring(s+1,e);
//				libraries[i]=libraries[i].substring(0,s);
//			}
//				
//		}
			
		
		
		return new ExperimentalDesign(libraries, samples, barcodes, labels);
	}

	private static double inferNumber(Pattern pat, String s) {
		Matcher m = pat.matcher(s);
		double re = 1;
		while (m.find()) {
			re = Double.parseDouble(m.group(1));
			// if directly followed by 4sU etc, return it immediately, otherwise use the last occurrence
			for (int i=0; i<regex_labels.length; i++) {
				Matcher m2 = regex_labels[i].matcher(s.substring(m.end()-1));
				if (m2.find() && m2.start()==0) return re;
			}
		}
		return re;
	}


	

	
	
	

}
