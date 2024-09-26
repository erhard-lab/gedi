package gedi.region.bam;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;
import java.util.Spliterator;
import java.util.Spliterators;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.Supplier;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gedi.bam.tools.BamUtils;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.ReferenceSequence;
import gedi.core.reference.Strand;
import gedi.core.reference.Strandness;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.sequence.SequenceProvider;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.FunctorUtils;
import gedi.util.FunctorUtils.MergeIterator;
import gedi.util.FunctorUtils.PeekIterator;
import gedi.util.StringUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.functions.FilteredSpliterator;
import gedi.util.functions.MappedSpliterator;
import gedi.util.io.text.LineWriter;
import gedi.util.mutable.MutableLong;
import gedi.util.sequence.DnaSequence;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;



/**
 * Due to softclipping in paired end mode with overlaps, the following may happen:
 * first read A and first read B come, then (maybe at a later pos) second read A
 * second read B may come with a greater pos (therefore A would be emitted here)
 * but if both A and B are overlappers (first and second read), they may still cover the same region!
 * 
 * Thus, the iterators may produce the same reference+genomicregion more than once!
 * If this is feed into a CIT, they will all be stored, so use CorrectCIT for that !
 * @author erhard
 *
 */
public class BamGenomicRegionStorage implements GenomicRegionStorage<AlignedReadsData> {
	
	
	public enum PairedEndHandling {
		AsMissingInformationIntron,JoinMates,IgnorePairedEnd,DropSecondRead,DropFirstRead
	}
	
	public enum PairedEndNoMateInRegionHandling {
		Drop,ReportSingle,Query
	}

	private Strandness[] strandness = null;
	private String[] fileNames;
	private SamReader[] files;
	private int[] cumNumCond; // number of conditions stored XR tag in files <=i

	private ContrastMapping mapping;

	private boolean ignoreVariations = false;
	private boolean removeMultiMapping = false;
	private boolean ignoreProperPair = false;
	
	private BiFunction<SAMRecord, SAMRecord, DnaSequence> barcode;
	private Predicate<SAMRecord> barcodeChecker;
	
	
//	private boolean joinMates = false;
//	private boolean ignorePairedEnd = false;
//	
//	private boolean forceIntersectionMates = false;
//	private boolean reportWithMissingMates = true;
	

	private PairedEndHandling pairedEndHandling = PairedEndHandling.AsMissingInformationIntron;
	private PairedEndNoMateInRegionHandling pairedEndNoMateInRegionHandling = PairedEndNoMateInRegionHandling.Drop;
	
	private boolean keepReadNames = false;
	private boolean onlyPrimary = false;
	
	private int minaqual = 0;

	private String name = "";
	
	private DynamicObject meta;
	
	private BamChecker bamChecker = null;
	private Predicate<SAMRecord> cbFilter;
	
	public BamGenomicRegionStorage(boolean strandspecific, BamMerge merge) {
		this(merge);
		setStrandness(strandspecific?Strandness.Sense:Strandness.Unspecific);
	}
	public BamGenomicRegionStorage(BamMerge merge) {
		this.fileNames = merge.getOriginalNames();
		name = merge.getName();

		
		this.files = new SamReader[fileNames.length];
		this.cumNumCond = new int[fileNames.length];
		Arrays.fill(cumNumCond, 1);
		for (int i=0; i<files.length; i++) {
			files[i]=obtainReader(i);
			for (String co : files[i].getFileHeader().getComments()) {
				if (co.startsWith("XR-count:"))
					throw new RuntimeException("XR not allowed when using merging!");
				else if (co.startsWith("@CO\tXR-count:"))
					throw new RuntimeException("XR not allowed when using merging!");
			}
		}
		ArrayUtils.cumSumInPlace(cumNumCond, 1);

		mapping = merge.getMapping();
		
		
		StringBuilder json = new StringBuilder().append("{\"conditions\" : [\n");
		for (int i=0; i<mapping.getNumMergedConditions(); i++) 
			json.append("\t{ \"name\" : \""+mapping.getMappedName(i)+"\"}").append(i+1<mapping.getNumMergedConditions()?",\n":"\n");
		meta = DynamicObject.parseJson(json.append("]}").toString());
		
		setStrandness(Strandness.Sense);
	}
	
	@Override
	public String getPath() {
		return fileNames[0];
	}

	public BamGenomicRegionStorage(boolean strandspecific, String... fileNames) throws IOException {
		this(fileNames);
		setStrandness(strandspecific?Strandness.Sense:Strandness.Unspecific);
	}
	public BamGenomicRegionStorage(String... fileNames) throws IOException {
		this.fileNames = fileNames;
		this.files = new SamReader[fileNames.length];
		this.cumNumCond = new int[fileNames.length];
		Arrays.fill(cumNumCond, 1);
		for (int i=0; i<files.length; i++) {
			files[i]=createReader(fileNames[i]);
			for (String co : files[i].getFileHeader().getComments()) {
				if (co.startsWith("XR-count:"))
					cumNumCond[i] = Integer.parseInt(co.substring("XR-count:".length()));
				else if (co.startsWith("@CO\tXR-count:"))
					cumNumCond[i] = Integer.parseInt(co.substring("@CO XR-count:".length()));
			}
		}
		ArrayUtils.cumSumInPlace(cumNumCond, 1);

		int lcs = fileNames.length>1?StringUtils.getLongestCommonSuffix(fileNames):0;

		// search for metadata
		ArrayList<DynamicObject> mc = new ArrayList<DynamicObject>();
		for (String f : fileNames) {
			File m = new File(f.replaceFirst(".bam$", ".metadata.json"));
			if (m.exists())
				try {
					mc.add(DynamicObject.parseJson(FileUtils.readAllText(m)));
				} catch (IOException e) {
					throw new RuntimeException("Could not read metadata file "+m,e);
				}
			else {
				StringBuilder json = new StringBuilder().append("{\"conditions\" : [ ")
							.append("\t{ \"name\" : \"")
							.append(getNameFromBam(f.substring(0,f.length()-lcs)))
							.append("\"}\n");
				mc.add(DynamicObject.parseJson(json.append("]}").toString()));
			}
		}
		meta = DynamicObject.merge(mc);
		
		StringBuilder sb = new StringBuilder();
		for (String f : fileNames){
			if (sb.length()>0) sb.append("_");
			sb.append(getNameFromBam(f));
		}
		name = sb.toString();
		
		if (new HashSet<>(Arrays.asList(getMetaDataConditions())).size()<getMetaDataConditions().length) {
			throw new RuntimeException("Did not read in bam file names properly: "+EI.wrap(getMetaDataConditions()).concat(",")+"(The logic here is: remove the common suffix from the paths, and, if 'Aligned' is part of the path, take only the part before it!)");
		}
		
		setStrandness(Strandness.Sense);
	}
	
	
	private String getNameFromBam(String n) throws IOException {
		if (new File(n).exists() && n.endsWith("possorted_genome_bam.bam")) // infer cell ranger name
			return new File(n).getCanonicalFile().getParentFile().getParentFile().getName();
		int idx = n.indexOf("Aligned");
		if (idx>0) n = n.substring(0,idx);
		n = StringUtils.removeFooter(n, ".bam");
		return n;
	}
//	private int[] readlengths = null;
//	public int[] getReadLengths() {
//		if (readlengths==null) {
//			readlengths = new int[isPairedEnd()?2:1];
//			SAMRecordIterator it = files[0].iterator();
//			while (it.hasNext() && ArrayUtils.min(readlengths)==0) {
//				SAMRecord rec = it.next();
//				int ri = !rec.getReadPairedFlag()||rec.getFirstOfPairFlag()?0:1;
//				if (readlengths[ri]==0) {
//					for (CigarElement ce : rec.getCigar().getCigarElements()) {
//						if (ce.getOperator().consumesReadBases())
//							readlengths[ri]+=ce.getLength();
//					}
//				}
//			}
//			it.close();
//		}
//		return readlengths;
//	}
	
	
	public void setBarcode(BiFunction<SAMRecord, SAMRecord, DnaSequence> barcode) {
		this.barcode = barcode;
	}
	
	public void setBarcode(String... attributes) {
		this.barcode = (r1, r2)->{
			StringBuilder sb = new StringBuilder();
			for (String a : attributes) {
				if (r1!=null) {
					String bc = r1.getStringAttribute(a);
					if (bc.indexOf('-')!=-1)
						bc = bc.substring(0, bc.indexOf('-'));
					sb.append(bc);
				}
					
				if (r2!=null){
					String bc = r2.getStringAttribute(a);
					if (bc.indexOf('-')!=-1)
						bc = bc.substring(0, bc.indexOf('-'));
					sb.append(bc);
				}
			}
			return new DnaSequence(sb);
		};
		
		this.barcodeChecker = (read)->{
			for (String a : attributes) {
				if (read!=null && read.getStringAttribute(a)==null)
					return false;
			}
			return true;
		};
	}
	
	public void setUmi() {
		this.barcode = (r1, r2)->{
			return new DnaSequence(StringUtils.removeFooter(StringUtils.splitField(r1.getReadName(), '#', 1),"_"));
		};
		
		this.barcodeChecker = (read)->{
			return read.getReadName().contains("#");
		};
		
		setMinimalAlignmentQuality(255);
	}
	
	public void setUmiOnt() {
		this.barcode = (r1, r2)->{
			String s = StringUtils.splitField(r1.getReadName(), '#', 0);
			int un = s.indexOf("_");
			return new DnaSequence(s.substring(0,un)+s.substring(un+1));
		};
		
		this.barcodeChecker = (read)->{
			int ha = read.getReadName().indexOf("#");
			int un = read.getReadName().indexOf("_");
			return un>=0 && ha>un;
		};
	}
	
	public void setUmi(String pattern) {
		Pattern pat = Pattern.compile(pattern);
		this.barcode = (r1, r2)->{
			Matcher m = pat.matcher(r1.getReadName());
			StringBuilder sb = new StringBuilder();
			for (int i=1; i<m.groupCount(); i++)
				sb.append(m.group(i));
			return new DnaSequence(sb.toString());
		};
		
		this.barcodeChecker = (read)->{
			return read.getReadName().contains("#");
		};
		
		setMinimalAlignmentQuality(255);
	}
	
	public void setDropseq(HashMap<String,String> filterBarcodes) {
		if (filterBarcodes!=null) 
			cbFilter = r->{
				String bc = r.getStringAttribute("XC");
				if (bc.indexOf('-')!=-1)
					bc = bc.substring(0, bc.indexOf('-'));
				return filterBarcodes.containsKey(bc);
			};
		
		setBarcode("XC","XM");
		setMinimalAlignmentQuality(255);
//		setIgnoreVariations(true);
	}
	
	public void set10x(HashMap<String,String> filterBarcodes) {
		if (filterBarcodes!=null) 
			cbFilter = r->{
				String bc = r.getStringAttribute("CB");
				if (bc.indexOf('-')!=-1)
					bc = bc.substring(0, bc.indexOf('-'));
				return filterBarcodes.containsKey(bc) && !r.getStringAttribute("UB").equals("-");
			};
		else {
			cbFilter = r->{
				String bc = r.getStringAttribute("CB");
				if (bc.indexOf('-')!=-1)
					bc = bc.substring(0, bc.indexOf('-'));
				return bc.length()>0 && !r.getStringAttribute("UB").equals("-");
			};
		}
		
		setBarcode("CB","UB");
		setMinimalAlignmentQuality(255);
//		setIgnoreVariations(true);
	}
	
	public void check(Genomic g, LineWriter inconsistent) {
		this.bamChecker = new BamChecker(g, inconsistent);
	}
	
	public String[] getFileNames() {
		return fileNames;
	}
	
	public void setMapping(ContrastMapping mapping) {
		this.mapping = mapping;
	}

	public ContrastMapping getMapping() {
		return mapping;
	}
	
	private SamReader createReader(String name) throws IOException {
		
		if (name.startsWith("http://") || name.startsWith("https://")) {
			String indexName = StringUtils.sha1(name+".bai")+".bai";
			if (!new File(indexName).exists()) {
				// download the index
				URL website = new URL(name+".bai");
				ReadableByteChannel rbc = Channels.newChannel(website.openStream());
				FileOutputStream fos = new FileOutputStream(indexName);
				fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
				fos.close();
			}
			
			return SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(SamInputResource.of(new URL(name)).index(new File(indexName)));
		}
		return SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(name));
	}

	@Override
	public Class<AlignedReadsData> getType() {
		return AlignedReadsData.class;
	}

	public String getName() {
		return name;
	}
	
	@Override
	public DynamicObject getMetaData() {
		return meta;
	}
	
	@Override
	public void setMetaData(DynamicObject meta) {
		this.meta = meta;
	}

	
	public BamGenomicRegionStorage setPairedEndHandling(PairedEndHandling pairedEndHandling) {
		this.pairedEndHandling = pairedEndHandling;
		return this;
	}
	
	public BamGenomicRegionStorage setPairedEndNoMateInRegionHandling(PairedEndNoMateInRegionHandling pairedEndNoMateInRegionHandling) {
		this.pairedEndNoMateInRegionHandling = pairedEndNoMateInRegionHandling;
		return this;
	}
	
	public BamGenomicRegionStorage detectVariations() {
		SAMRecordIterator it = files[0].iterator();
		SAMRecord r = it.next();
		it.close();
		return setIgnoreVariations(r.getStringAttribute("XV")==null && r.getStringAttribute("MD")==null);
	}
	
	public BamGenomicRegionStorage setIgnoreVariations(boolean ignore) {
		this.ignoreVariations = ignore;
		return this;
	}
	public BamGenomicRegionStorage setOnlyPrimary(boolean onlyPrimary) {
		this.onlyPrimary = onlyPrimary;
		return this;
	}
	
	public BamGenomicRegionStorage setKeepReadNames(boolean keepReadNames) {
		this.keepReadNames = keepReadNames;
		return this;
	}

	public BamGenomicRegionStorage setReportWithMissingMates(boolean reportWithMissingMates) {
//		this.reportWithMissingMates = reportWithMissingMates;
		this.pairedEndNoMateInRegionHandling = reportWithMissingMates?PairedEndNoMateInRegionHandling.ReportSingle:PairedEndNoMateInRegionHandling.Drop;
		return this;
	}
	
	public BamGenomicRegionStorage setForceIntersectionMates(boolean forceIntersectionMates) {
//		this.forceIntersectionMates = forceIntersectionMates;
		this.pairedEndNoMateInRegionHandling = forceIntersectionMates?PairedEndNoMateInRegionHandling.Query:PairedEndNoMateInRegionHandling.Drop;

		return this;
	}
	
	public BamGenomicRegionStorage setRemoveMultiMapping(boolean remove) {
		this.removeMultiMapping = remove;
		return this;
	}

	public BamGenomicRegionStorage setMinimalAlignmentQuality(int minaqual) {
		this.minaqual  = minaqual;
		return this;
	}
	
	public int getMinimalAlignmentQuality() {
		return minaqual;
	}
	
	public BamGenomicRegionStorage setIgnoreProperPair(boolean ignoreProperPair) {
		this.ignoreProperPair = ignoreProperPair;
		return this;
	}
	
	
	public BamGenomicRegionStorage setJoinMates(boolean join) {
//		this.joinMates = join;
		this.pairedEndHandling = join?PairedEndHandling.JoinMates:PairedEndHandling.AsMissingInformationIntron;
		return this;
	}

	public BamGenomicRegionStorage setIgnorePairedEnd(boolean ignoreMates) {
//		this.ignorePairedEnd = ignoreMates;
		this.pairedEndHandling = ignoreMates?PairedEndHandling.IgnorePairedEnd:PairedEndHandling.AsMissingInformationIntron;
		return this;
	}
	
	private enum PairedEnd {
		Unknown,Paired,Unpaired;
	}
	private PairedEnd pairedEndInfo = PairedEnd.Unknown;
	public boolean isPairedEnd() {
		if (pairedEndInfo==PairedEnd.Unknown) {
			SAMRecordIterator it = files[0].iterator();
			pairedEndInfo = it.next().getReadPairedFlag()?PairedEnd.Paired:PairedEnd.Unpaired;
			it.close();
		}
		return pairedEndInfo==PairedEnd.Paired;
	}
	
	public boolean isStrandSpecific(int file) {
		return strandness[file]!=Strandness.Unspecific;
	}
	
	public boolean isAnyStrandSpecific() {
		for (int i=0; i<files.length; i++)
			if (isStrandSpecific(i))
				return true;
		return false;
	}
	
	public BamGenomicRegionStorage setStrandness(Strandness strandness) {
		this.strandness = new Strandness[files.length];
		for (int i=0; i<files.length; i++)
			setStrandness(i,strandness);
		return this;
	}
	
	public BamGenomicRegionStorage setStrandness(Strandness[] strandness) {
		for (int i=0; i<fileNames.length; i++)
			setStrandness(i,strandness[i]);
		return this;
	}
	
	public BamGenomicRegionStorage setStrandness(int file, Strandness strandness) {
		if (strandness.equals(Strandness.AutoDetect))
			throw new RuntimeException("Not implemented!");
		this.strandness[file] = strandness;
		return this;
	}
	public Strandness getStrandness(int file) {
		return strandness[file];
	}

	
	private HashMap<ReferenceSequence,Integer> refs=null;
	private ReentrantLock refLock = new ReentrantLock();
	@Override
	public Set<ReferenceSequence> getReferenceSequences() {
		if (refs==null) {
			refLock.lock();
			if (refs==null) {
				boolean spec = false;
				boolean unspec = false;
				for (int i=0; i<files.length; i++)
					if (isStrandSpecific(i))
						spec = true;
					else
						unspec = true;
				
				HashMap<ReferenceSequence,Integer> refs = new HashMap<>();
				for (SamReader s : files)
					for (SAMSequenceRecord seq : s.getFileHeader().getSequenceDictionary().getSequences()) {
						if (spec) {
							refs.put(Chromosome.obtain(seq.getSequenceName(),Strand.Plus),seq.getSequenceLength());
							refs.put(Chromosome.obtain(seq.getSequenceName(),Strand.Minus),seq.getSequenceLength());
						} 
						if (unspec)
							refs.put(Chromosome.obtain(seq.getSequenceName(),Strand.Independent),seq.getSequenceLength());
					}
				this.refs = refs;
			}
			refLock.unlock();
		}
		
		return refs.keySet();
	}

	private int getSequenceLength(ReferenceSequence ref) {
		getReferenceSequences();
		Integer re = refs.get(ref);
		return re==null?0:re;
	}


	public void close() {
		for (SamReader s : files)
			try {
				s.close();
			} catch (IOException e) {
			}
	}


	public Spliterator<GenomicRegion> iterateGenomicRegions(ReferenceSequence ref) {
		return new BufferingGenomicRegionSpliterator(ref);
//		return convert?new ConvertedGenomicRegionSpliterator(ref):new GenomicRegionSpliterator(ref);
	}

	public Spliterator<GenomicRegion> iterateGenomicRegions(ReferenceSequence ref, int start, int end) {
		return new BufferingGenomicRegionSpliterator(ref, start, end);
//		return convert?new ConvertedGenomicRegionSpliterator(ref, start, end):new GenomicRegionSpliterator(ref, start, end);
	}

	public AlignedReadsData getAlignedReadsData(ReferenceSequence ref, GenomicRegion region) {
		return getAlignedReadsData(ref, region, null);
	}
	public AlignedReadsData getAlignedReadsData(ReferenceSequence ref, GenomicRegion region, SequenceProvider sequenceProvider) {

		if (region instanceof FactoryGenomicRegion) {
			AlignedReadsData re2 = ((FactoryGenomicRegion)region).create();
			if (mapping==null) return re2;
			return new ConditionMappedAlignedReadsData(re2,mapping);
		}
			
		

//		if (region instanceof ConditionMappedAlignedReadsData && ((ConditionMappedAlignedReadsData)region).getParent() instanceof FactoryGenomicRegion)
//			return ((FactoryGenomicRegion)((ConditionMappedAlignedReadsData)region).getParent()).create();

		ImmutableReferenceGenomicRegion<AlignedReadsData> re = ei(new ImmutableReferenceGenomicRegion<>(ref, region)).filter(r->region.equals(r.getRegion())).getUniqueResult(true, false);
		return re==null?null:re.getData();
		

//		SAMRecord[][] re = getRecords(ref, region);
//		if (re==null) return null;
//
//		CharSequence seq = null;
//		if (sequenceProvider!=null)
//			seq = sequenceProvider.getPlusSequence(ref.getName(), region);
//		
//		BamAlignedReadsData re2 = new BamAlignedReadsData(re,cumNumCond,ref.getStrand()==Strand.Minus,seq);
//
//		if (mapping==null) return re2;
//
//		return new ConditionMappedAlignedReadsData(re2, mapping);
	}


//	/**
//	 * Returns null, if this region is not present in any sam file!
//	 * @param ref
//	 * @param region
//	 * @return
//	 */
//	public SAMRecord[][] getRecords(ReferenceSequence ref, GenomicRegion region) {
//		if (region instanceof RecordsGenomicRegion) 
//			return ((RecordsGenomicRegion)region).getRecords();
//		if (region instanceof ConditionMappedAlignedReadsData && ((ConditionMappedAlignedReadsData)region).getParent() instanceof RecordsGenomicRegion)
//			((RecordsGenomicRegion)((ConditionMappedAlignedReadsData)region).getParent()).getRecords();
//
//		SAMRecord[][] re = null;
//		ArrayList<SAMRecord> c = new ArrayList<SAMRecord>();
//		for (int i=0; i<files.length; i++) {
//			SAMRecordIterator it = files[i].query(ref.getName(), region.getStart()+1,region.getEnd(),false);
//			while (it.hasNext()) {
//				SAMRecord n = it.next();
//				if ((!strandSpecific || BamUtils.isValidStrand(ref.getStrand(), n)) && BamUtils.getRecordsGenomicRegion(n).equals(region))
//					c.add(n);
//			}
//			it.close();
//
//			if (!c.isEmpty()) {
//				if (re==null)
//					re = new SAMRecord[files.length][];
//				re[i] = c.toArray(new SAMRecord[0]);
//				c.clear();
//			}
//		}
//		if (re!=null)
//			for (int i=0; i<re.length; i++)
//				if (re[i]==null) 
//					re[i] = new SAMRecord[0];
//		return re;
//	}
	
//	LinkedBlockingQueue<SamReader>[] readers;
	
	private SamReader obtainReader(int file) {
//		if (readers==null) readers = new LinkedBlockingQueue[fileNames.length];
//		if (readers[file]==null) readers[file] = new LinkedBlockingQueue<SamReader>();
//		
//		LinkedBlockingQueue<SamReader> ll = readers[file];
//
//		SamReader re = ll.poll();
		SamReader re = null;
//		for (;!ll.isEmpty(); re = ll.poll())
//			if (!ll.isEmpty()) {
//				re.close();
//				System.out.println("closed reader for "+fileNames[file]);
//			}
		
		if (re==null) {
//			System.out.println("Creating new reader for "+fileNames[file]);
			try {
				re = createReader(fileNames[file]);
			} catch (IOException e) {
				throw new RuntimeException("Cannot read bam file!",e);
			}
//			if (false) {
//				try {
//					Object reader = ReflectionUtils.get(re, "mReader");
//					MemoryBAMFileIndex mIndex = new MemoryBAMFileIndex(ReflectionUtils.get(reader, "mIndexFile"), re.getFileHeader().getSequenceDictionary(), true);
//					ReflectionUtils.set(reader, "mIndex", mIndex);
//				} catch (Exception e) {
//					throw new RuntimeException("Cannot load index!",e);
//				}
//			}
//			re.enableIndexCaching(true);
		}
		
//		System.out.println("Reusing reader for "+fileNames[file]);
			
		return re;
	}
	
	private void returnReader(int file, SamReader sf) {
//		if (readers==null) readers = new LinkedBlockingQueue[fileNames.length];
//		if (readers[file]==null) readers[file] = new LinkedBlockingQueue<SamReader>();
//		try {
////			System.out.println("Returning reader for "+fileNames[file]);
//			readers[file].put(sf);
//		} catch (InterruptedException e) {
//			throw new RuntimeException(e);
//		}
		try {
			sf.close();
		} catch (IOException e) {
		}
//		sf.close();
	}

	public boolean contains(ReferenceSequence ref, GenomicRegion region) {
		for (int i=0; i<files.length; i++) {
			SAMRecordIterator it = files[i].query(ref.getName(), region.getStart()+1,region.getEnd(),false);
			while (it.hasNext()) {
				SAMRecord n = it.next();
				if (BamUtils.isValidStrand(ref.getStrand(),strandness[i], n, null) && BamUtils.getRecordsGenomicRegion(n).equals(region)) {
					it.close();
					return true;
				}
			}
			it.close();
		}
		return false;
	}

//	
//	@Deprecated
//	private class GenomicRegionSpliterator implements Spliterator<GenomicRegion> {
//
//		private static final int START_ESTIMATE = 10;
//		private int minSpanForSplit=10000;
//
//		private ReferenceSequence ref;
//		private int pos;
//		private int end;
//
//		private int regionsReturned = 0;
//		private int basesCovered = 0;
//
//		private SAMFileReader[] files;
//		private Iterator<RecordsGenomicRegion> iterator;
//
//		public GenomicRegionSpliterator(ReferenceSequence ref) {
//			this(ref,0,getSequenceInfo(ref).getSequenceLength());
//		}
//
//		public GenomicRegionSpliterator(ReferenceSequence ref, int start, int end) {
//			this.ref = ref;
//			pos = start;
//			this.end = end;
//			files = new SAMFileReader[BamGenomicRegionStorage.this.files.length];
//
//			Iterator<RecordsGenomicRegion>[] its = new Iterator[files.length];
//			for (int i=0; i<its.length; i++) {
//				files[i] = obtainReader(i);
//				files[i].setValidationStringency(ValidationStringency.SILENT);
//				//				files[i]= BamGenomicRegionStorage.this.files[i];
//				SAMRecordIterator raw = files[i].query(ref.getName(), pos+1, end, false);
//				Iterator<SAMRecord> filtered = strandSpecific?FunctorUtils.filteredIterator(raw,a -> BamUtils.isValidStrand(ref.getStrand(),a)):raw;
//				MappedIterator<SAMRecord,RecordsGenomicRegion> mapped = FunctorUtils.mappedIterator(filtered, BamUtils::getRecordsGenomicRegion);
//				Iterator<RecordsGenomicRegion> resort = FunctorUtils.resortIterator(mapped, new GenomicRegion.StartPositionComparator(), getComparator());
//				Iterator<RecordsGenomicRegion> unified = FunctorUtils.reduceIterator(resort, getComparator(), (a,b)->a.merge(b));
//
//				its[i] = unified;
//			}
//			iterator = FunctorUtils.mappedIterator(FunctorUtils.parallellIterator(its, getComparator(), RecordsGenomicRegion.class),(a)->{
//				int i=0;
//				for (;i<a.length && a[i]==null; i++);
//				if (i==a.length) return null;
//				return a[i].join(a);
//			});
//			//FunctorUtils.unifyIterator(FunctorUtils.mergeIterator(its, getComparator()),getComparator());
//		}
//
//		@Override
//		public Comparator<? super GenomicRegion> getComparator() {
//			return (GenomicRegion a, GenomicRegion b) -> a.compareTo(b);
//		}
//
//
//		@Override
//		public boolean tryAdvance(Consumer<? super GenomicRegion> action) {
//			if (iterator.hasNext()) {
//				GenomicRegion re = iterator.next();
//
//				basesCovered+=re.getStart()-pos;
//				regionsReturned++;
//
//				pos = re.getStart();
//				if (re.getStart()>=end) return false;
//				action.accept(re);
//				return true;
//			}
//			for (int i=0; i<files.length; i++)
//				files[i].close();
//			return false;
//		}
//
//		@Override
//		public Spliterator<GenomicRegion> trySplit() {
//			if (!iterator.hasNext()) return null;
//			if (end-pos<minSpanForSplit) return null;
//
//			int center = (pos+end)/2;
//			if (center==pos) return null;
//
//			GenomicRegionSpliterator re = new GenomicRegionSpliterator(ref,center,end);
//			re.regionsReturned = regionsReturned;
//			re.basesCovered = basesCovered;
//			end = center;
//			return re;			
//		}
//
//		@Override
//		public long estimateSize() {
//			if (regionsReturned<START_ESTIMATE) return Long.MAX_VALUE;
//			return (long) (Math.abs(end-pos)*(regionsReturned/(double)basesCovered));
//		}
//
//		@Override
//		public int characteristics() {
//			return DISTINCT|NONNULL|SORTED|ORDERED|IMMUTABLE;
//		}
//
//	}

	


//	private class ConvertedGenomicRegionSpliterator implements Spliterator<GenomicRegion> {
//
//		private static final int START_ESTIMATE = 10;
//		private int minSpanForSplit=10000;
//
//		private ReferenceSequence ref;
//		private int pos;
//		private int end;
//
//		private int regionsReturned = 0;
//		private int basesCovered = 0;
//
//		private SAMFileReader[] files;
//		private SAMFileReader[] mateFiles;
//		private PeekIterator<SAMRecord> defIterator;
//		private SAMRecordIterator[] raw;
//		
//		private LinkedList<FactoryGenomicRegion> queue = new LinkedList<FactoryGenomicRegion>();
//
//
//		public ConvertedGenomicRegionSpliterator(ReferenceSequence ref) {
//			this(ref,0,getSequenceInfo(ref).getSequenceLength());
//		}
//
//		public ConvertedGenomicRegionSpliterator(ReferenceSequence ref, int start, int end) {
//			this.ref = ref;
//			pos = start;
//			this.end = end;
//			files = new SAMFileReader[BamGenomicRegionStorage.this.files.length];
//			raw = new SAMRecordIterator[BamGenomicRegionStorage.this.files.length];
//			
//			if (isPairedEnd() && !ignorePairedEnd)
//				cacheMates(ref.getName(), start, end);
//		}
//
//		private HashMap<String,SAMRecord>[] mateCache;
//		private void cacheMates(String name, int start, int end) {
//			mateFiles = new SAMFileReader[BamGenomicRegionStorage.this.files.length];
//			mateCache = new HashMap[files.length];
//			for (int i=0; i<files.length; i++) {
//				mateFiles[i] = obtainReader(i);
//				mateFiles[i].setValidationStringency(ValidationStringency.SILENT);
//				mateCache[i] = new HashMap<String, SAMRecord>();
//				try {
//					SAMRecordIterator it = mateFiles[i].query(name, start, end,false);
//	//				System.out.println(fileNames[i]+"\t"+name+":"+start+"-"+end);
//					try {
//						while (it.hasNext()) {
//							final SAMRecord next = it.next();
//							if (next.getSecondOfPairFlag()) {
//								mateCache[i].put(BamUtils.getPairId(next), next);
//							}
//						}
//	//					System.out.println(mateCache[i].size());
//					} finally {
//						it.close();
//					}
//				} catch (NullPointerException e) {
//					mateFiles[i].query(name, start, end,false);
//				}
//			}
//		}
//
//
//		private SAMRecord queryMate(int file, SAMRecord rec, boolean onlyCached) {
//			
//			if (mateCache!=null) {
//				SAMRecord mate = mateCache[file].get(BamUtils.getPairId(rec));
//				if (mate!=null && mate.getReferenceName().equals(rec.getMateReferenceName()) && mate.getAlignmentStart()==rec.getMateAlignmentStart()) {
//					if (minaqual>0 && mate.getMappingQuality()<minaqual) return null;
//					if (removeMultiMapping && mate.getIntegerAttribute("NH")!=null && mate.getIntegerAttribute("NH")>1) return null;
//					return mate;
//				}
//			}
//			if (onlyCached)
//				return null;
//			
//			SAMRecordIterator it = mateFiles[file].queryAlignmentStart(rec.getMateReferenceName(), rec.getMateAlignmentStart());
//			try {
//				while (it.hasNext()) {
//					final SAMRecord next = it.next();
//					if (next.getFirstOfPairFlag()) continue;
//					if (BamUtils.getPairId(rec).equals(BamUtils.getPairId(next))) {
//						if ((minaqual==0 || next.getMappingQuality()>=minaqual) && !(removeMultiMapping && next.getIntegerAttribute("NH")!=null && next.getIntegerAttribute("NH")>1)) 
//							return next;
//					}
//				}
//				return null;
//			} finally {
//				it.close();
//			}
//		}
//		
//		private PeekIterator<SAMRecord> iterator() {
//			if (defIterator==null) {
//				Iterator<SAMRecord>[] its = new Iterator[files.length];
//				for (int i=0; i<its.length; i++) {
//					files[i] = obtainReader(i);
//					files[i].setValidationStringency(ValidationStringency.SILENT);
//					
//					//				files[i]= BamGenomicRegionStorage.this.files[i];
//					raw[i] = files[i].query(ref.getName(), pos+1, end, false);
//					Iterator<SAMRecord> filtered = strandSpecific?FunctorUtils.filteredIterator(raw[i],
//							a -> BamUtils.isValidStrand(ref.getStrand(),a)):raw[i];
//
//					if (minaqual>0)
//						filtered = FunctorUtils.filteredIterator(filtered,
//							a -> a.getMappingQuality()>=minaqual);
//					if (removeMultiMapping)
//						filtered = FunctorUtils.filteredIterator(filtered,
//								a -> a.getIntegerAttribute("NH")==null || a.getIntegerAttribute("NH")==1);
//					
//
//					// filter first mate pair here if necessary
//					if (isPairedEnd() && !ignorePairedEnd)
//						filtered = FunctorUtils.filteredIterator(filtered,
//								a -> a.getFirstOfPairFlag());
//
//					its[i] = filtered;
//				}
//				defIterator = FunctorUtils.peekIterator( 
//						FunctorUtils.mergeIterator(its, (a,b)->{
//							return Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
//						}
//								));
//			}
//			return defIterator;
//		}
//
//		@Override
//		public Comparator<? super GenomicRegion> getComparator() {
//			return (GenomicRegion a, GenomicRegion b) -> a.compareTo(b);
//		}
//
//
//		@Override
//		public boolean tryAdvance(Consumer<? super GenomicRegion> action) {
//			if (!queue.isEmpty()) {
//				action.accept(queue.removeFirst());
//				return true;
//			}
//			PeekIterator<SAMRecord> iterator = iterator();
//			advance:while (iterator.hasNext()) {
//				SAMRecord rec = iterator.next();
//				
//				basesCovered+=rec.getAlignmentStart()-1-pos;
//				regionsReturned++;
//
//				pos = rec.getAlignmentStart()-1;
//				if (rec.getAlignmentStart()-1>=end) return false;
//
//				HashMap<FactoryGenomicRegion,FactoryGenomicRegion> map = new HashMap<FactoryGenomicRegion, FactoryGenomicRegion>();
//				do {
//					if (!map.isEmpty()) rec = iterator.next();
//
////					if (rec.getReadName().equals("3845271"))
////						System.out.println();
//
//					if (rec.getReadPairedFlag() && !ignorePairedEnd ) {
//						// paired end is super inefficient (stupid samtools!)
//
//						// this is already done by the iterator!
//						//if (rec.getFirstOfPairFlag()) {
//
//						int file = ((MergeIterator<SAMRecord>)iterator.getParent()).getIteratorIndex();
//						SAMRecord mate = queryMate(file,rec,true);
//						if (mate==null || (!ignoreProperPair && !mate.getProperPairFlag()))
//							continue;
//						FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,mate,cumNumCond, joinMates,false,ignoreVariations,false);
//						if (!fac.isConsistent()) continue;
//
//						if (!map.containsKey(fac)) map.put(fac,fac);
//						else fac = map.get(fac);
//
//						fac.add(rec, mate, file);
//						//}
//
//					} else {
//
//						FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,cumNumCond,ignoreVariations,false);
//						if (!map.containsKey(fac)) map.put(fac,fac);
//						else fac = map.get(fac);
//
//						int file = ((MergeIterator<SAMRecord>)iterator.getParent()).getIteratorIndex();
//						fac.add(rec, file);
//					}
//
//				} while (iterator.hasNext() && iterator.peek().getAlignmentStart()-1==pos && !map.isEmpty());
//
//				if (!map.isEmpty()) {
//					FactoryGenomicRegion[] arr = map.keySet().toArray(new FactoryGenomicRegion[0]);
//					Arrays.sort(arr,getComparator());
//					for (int i=0; i<arr.length; i++)
//						queue.add(arr[i]);
//	
//					action.accept(queue.removeFirst());
//					return true;
//				}
//			}
//			for (int i=0; i<files.length; i++) {
//				returnReader(i,files[i]);
//				raw[i].close();
////				files[i].close();
//				if (mateFiles!=null)
//					returnReader(i,mateFiles[i]);
////					mateFiles[i].close();
//			}
//			return false;
//		}
//
//
//
//		@Override
//		public Spliterator<GenomicRegion> trySplit() {
////			if (!iterator().hasNext()) return null;
////			if (end-pos<minSpanForSplit) return null;
////
////			int center = (pos+end)/2;
////			if (center==pos) return null;
////
////			GenomicRegionSpliterator re = new GenomicRegionSpliterator(ref,center,end,false);
////			re.regionsReturned = regionsReturned;
////			re.basesCovered = basesCovered;
////			end = center;
////			return re;
//			return null;
//		}
//
//		@Override
//		public long estimateSize() {
//			if (regionsReturned<START_ESTIMATE) return Long.MAX_VALUE;
//			return (long) (Math.abs(end-pos)*(regionsReturned/(double)basesCovered));
//		}
//
//		@Override
//		public int characteristics() {
//			return DISTINCT|NONNULL|SORTED|ORDERED|IMMUTABLE;
//		}
//
//	}
	
	private static class SAMRecordList {
		SAMRecord r;
		int file;
		SAMRecordList next;
		public SAMRecordList(SAMRecord r, SAMRecordList next, int file) {
			super();
			this.r = r;
			this.next = next;
			this.file = file;
		}
		
	}
	
	private class BufferingGenomicRegionSpliterator implements Spliterator<GenomicRegion> {

		private static final int START_ESTIMATE = 10;

		private ReferenceSequence ref;
		private int pos;
		private int end;

		private int regionsReturned = 0;
		private int basesCovered = 0;

		private SamReader[] files;
		private PeekIterator<SAMRecord> defIterator;
		private SAMRecordIterator[] raw;
		
		private LinkedList<FactoryGenomicRegion> queue = new LinkedList<FactoryGenomicRegion>();


		public BufferingGenomicRegionSpliterator(ReferenceSequence ref) {
			this(ref,0,getSequenceLength(ref));
		}

		public BufferingGenomicRegionSpliterator(ReferenceSequence ref, int start, int end) {
			this.ref = ref;
			pos = start;
			this.end = end;
			files = new SamReader[BamGenomicRegionStorage.this.files.length];
			raw = new SAMRecordIterator[BamGenomicRegionStorage.this.files.length];
		}


		
		private PeekIterator<SAMRecord> iterator(int start, int end) {
			if (defIterator==null) {
				Iterator<SAMRecord>[] its = new Iterator[files.length];
				for (int i=0; i<its.length; i++) {
					files[i] = obtainReader(i);
					
					//				files[i]= BamGenomicRegionStorage.this.files[i];
					String name = ref.getName();
					if (files[i].getFileHeader().getSequenceIndex(name)==-1 && files[i].getFileHeader().getSequenceIndex(StringUtils.removeHeader(name,"chr"))!=-1)
						name = StringUtils.removeHeader(name,"chr");
					if (files[i].getFileHeader().getSequenceIndex(name)==-1 && files[i].getFileHeader().getSequenceIndex("chr"+name)!=-1)
						name = "chr"+name;

					if (files[i].getFileHeader().getSequenceIndex(name)==-1 && name.equals("MT") && files[i].getFileHeader().getSequenceIndex("M")!=-1)
						name = "M";
					if (files[i].getFileHeader().getSequenceIndex(name)==-1 && name.equals("MT") && files[i].getFileHeader().getSequenceIndex("chrM")!=-1)
						name = "chrM";
					
					if (files[i].getFileHeader().getSequenceIndex(name)==-1) {
						raw[i] = null;
						its[i] = EI.empty();
					}
					else {
						raw[i] = files[i].query(name, start, end, false);
					
					
						Iterator<SAMRecord> filtered = raw[i];
	//							strandSpecific?FunctorUtils.filteredIterator(raw[i],
	//							a -> BamUtils.isValidStrand(ref.getStrand(),a)):raw[i];
						// for pairedEnd: handle that differently!
						
						if (barcodeChecker!=null)
							filtered = FunctorUtils.filteredIterator(filtered, barcodeChecker);
						
						if (minaqual>0)
							filtered = FunctorUtils.filteredIterator(filtered,
								a -> a.getMappingQuality()>=minaqual);
						if (removeMultiMapping)
							filtered = FunctorUtils.filteredIterator(filtered,
									a -> a.getIntegerAttribute("NH")==null || a.getIntegerAttribute("NH")==1);
						
						if (cbFilter!=null)
							filtered = FunctorUtils.filteredIterator(filtered, cbFilter);
	
						filtered = FunctorUtils.filteredIterator(filtered,
								a -> !a.getReadUnmappedFlag());
	
					
						its[i] = filtered;
					}
				}
				defIterator = FunctorUtils.peekIterator( 
						FunctorUtils.mergeIterator(its, (a,b)->{
							return Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
						}
								));
			}
			return defIterator;
		}

		@Override
		public Comparator<? super GenomicRegion> getComparator() {
			return (GenomicRegion a, GenomicRegion b) -> a.compareTo(b);
		}

		private HashMap<String,SAMRecordList> mateBuffer = new HashMap<String, SAMRecordList>();
		private boolean done = false;
		@Override
		public boolean tryAdvance(Consumer<? super GenomicRegion> action) {
			try {
				if (!queue.isEmpty()) {
					action.accept(queue.removeFirst());
					return true;
				}
				if (!done) {
					PeekIterator<SAMRecord> iterator = iterator(pos+1, end);
					while (iterator.hasNext()) {
						SAMRecord rec = iterator.peek();
						
						basesCovered+=rec.getAlignmentStart()-1-pos;
						regionsReturned++;
		
						pos = rec.getAlignmentStart()-1;
						if (rec.getAlignmentStart()-1>=end) return false;
		
						HashMap<FactoryGenomicRegion,FactoryGenomicRegion> map = new HashMap<FactoryGenomicRegion, FactoryGenomicRegion>();
						do {
							rec = iterator.next();
							processRecord(true,rec, iterator, map);
		
						} while (iterator.hasNext() && iterator.peek().getAlignmentStart()-1==pos && !map.isEmpty());
		
						if (!map.isEmpty()) {
							FactoryGenomicRegion[] arr = map.keySet().toArray(new FactoryGenomicRegion[0]);
							Arrays.sort(arr,getComparator());
							for (int i=0; i<arr.length; i++)
								queue.add(arr[i]);
			
							action.accept(queue.removeFirst());
							return true;
						}
					}
					done = true;
				}
				for (int i=0; i<files.length; i++) {
					returnReader(i,files[i]);
					if (raw[i]!=null)
						raw[i].close();
	//				files[i].close();
				}
			} catch (htsjdk.samtools.util.RuntimeIOException e) {
				return false;
			}
			
//			for (SAMRecordList l : mateBuffer.values())
//				for (; l!=null; l=l.next)
//					System.out.print("From Matebuffer: "+l.r.getSAMString());
			
			// query mates that are not within pos-end
			int min = -1;
			int max = -1;
//			if (forceIntersectionMates) {
			if (pairedEndNoMateInRegionHandling==PairedEndNoMateInRegionHandling.Query) {
				for (SAMRecordList l : mateBuffer.values())
					for (; l!=null; l=l.next) {
						if (min==-1) {
							max = min = l.r.getMateAlignmentStart();
						} else {
							min = Math.min(l.r.getMateAlignmentStart(), min);
							max = Math.max(l.r.getMateAlignmentStart(), max);
						}
					}
			}
			if (min>-1) {
				defIterator = null;
				PeekIterator<SAMRecord> iterator = iterator(min, max);
				while (iterator.hasNext()) {
					SAMRecord rec = iterator.peek();
					
					basesCovered+=rec.getAlignmentStart()-1-pos;
					regionsReturned++;
	
					pos = rec.getAlignmentStart()-1;
					
					HashMap<FactoryGenomicRegion,FactoryGenomicRegion> map = new HashMap<FactoryGenomicRegion, FactoryGenomicRegion>();
					do {
						rec = iterator.next();
						processRecord(false, rec, iterator, map);
	
					} while (iterator.hasNext() && iterator.peek().getAlignmentStart()-1==pos && !map.isEmpty());
	
					if (!map.isEmpty()) {
						FactoryGenomicRegion[] arr = map.keySet().toArray(new FactoryGenomicRegion[0]);
						Arrays.sort(arr,getComparator());
						for (int i=0; i<arr.length; i++)
							queue.add(arr[i]);
					}
					map.clear();
				}
				for (int i=0; i<files.length; i++) {
					returnReader(i,files[i]);
					if (raw[i]!=null)
						raw[i].close();
//					files[i].close();
				}

				mateBuffer.clear();
				if (!queue.isEmpty()) {
					action.accept(queue.removeFirst());
					return true;
				}
			}
			
//			if (reportWithMissingMates) {
			if (pairedEndNoMateInRegionHandling==PairedEndNoMateInRegionHandling.ReportSingle) {	
				HashMap<FactoryGenomicRegion,FactoryGenomicRegion> map = new HashMap<FactoryGenomicRegion, FactoryGenomicRegion>();
				for (SAMRecordList ll : mateBuffer.values()) {
					for (; ll!=null; ll=ll.next){
						
						if (isAnyStrandSpecific()){
							SAMRecord first = ll.r.getFirstOfPairFlag()?ll.r:null;
							SAMRecord second = ll.r.getFirstOfPairFlag()?null:ll.r;
							if (BamUtils.isValidStrand(ref.getStrand(), strandness[ll.file], first,second)) {
								FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(ll.r,null,cumNumCond, pairedEndHandling==PairedEndHandling.JoinMates,false,ignoreVariations,keepReadNames,barcode);
								if (!map.containsKey(fac)) map.put(fac,fac);
								else fac = map.get(fac);
								
								fac.add(ll.r, ll.file);

							}
						}
						
					}
				}
				for (FactoryGenomicRegion r : map.keySet())
					queue.add(r);
				mateBuffer.clear();
				
				if (!queue.isEmpty()) {
					action.accept(queue.removeFirst());
					return true;
				}
			}
			
//			System.out.println(ref);
//			System.out.println(EI.wrap(statslab).concat(","));
//			System.out.println(EI.wrap(stats).concat(","));
			
			return false;
		}



		private void processRecord(boolean fillMateBuffer, SAMRecord rec, PeekIterator<SAMRecord> iterator, HashMap<FactoryGenomicRegion, FactoryGenomicRegion> map) {
//			System.out.print(rec.getSAMString());
			if (rec.getNotPrimaryAlignmentFlag() && onlyPrimary)
				return;
			
			if (rec.getReadPairedFlag() && pairedEndHandling==PairedEndHandling.DropSecondRead) {
				if (rec.getFirstOfPairFlag())
					processSingleEnd(fillMateBuffer, rec, iterator, map);
				
			} 
			else if (rec.getReadPairedFlag() && pairedEndHandling==PairedEndHandling.DropFirstRead) {
				if (!rec.getFirstOfPairFlag())
					processSingleEnd(fillMateBuffer, rec, iterator, map);
				
			} 
			else if (rec.getReadPairedFlag() && pairedEndHandling!=PairedEndHandling.IgnorePairedEnd) {
				processPairedEnd(fillMateBuffer, rec, iterator, map);

			} else {
				processSingleEnd(fillMateBuffer, rec, iterator, map);
			}

		}
		String[] statslab= {"SE incoming", "SE added", "PE incoming (read)", "PE added (read)", "PE inconsistent (pair)", "PE consistent (pair)"};
		int[] stats = new int[2+4];
		private void processSingleEnd(boolean fillMateBuffer, SAMRecord rec, PeekIterator<SAMRecord> iterator, HashMap<FactoryGenomicRegion, FactoryGenomicRegion> map) {
			int file = ((MergeIterator<SAMRecord>)iterator.getParent()).getIteratorIndex();
			stats[0]++;
			
			if (!BamUtils.isValidStrand(ref.getStrand(), strandness[file], rec,null))
				return;
			
			if (!ignoreVariations && bamChecker!=null) {
				FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,cumNumCond,ignoreVariations,keepReadNames, barcode);
				fac.add(rec, file);
				bamChecker.check(rec,new ImmutableReferenceGenomicRegion<>(ref, fac,fac.create()));
			}

			FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,cumNumCond,ignoreVariations,keepReadNames,barcode);
			if (!map.containsKey(fac)) map.put(fac,fac);
			else fac = map.get(fac);
			
			stats[1]++;
			fac.add(rec, file);
		}
		
		private void processPairedEnd(boolean fillMateBuffer, SAMRecord rec, PeekIterator<SAMRecord> iterator, HashMap<FactoryGenomicRegion, FactoryGenomicRegion> map) {
			stats[2]++;
			
			// paired end is super inefficient (stupid samtools!)
			if (!ignoreProperPair && !rec.getProperPairFlag()) return;

			stats[3]++;
			
			int file = ((MergeIterator<SAMRecord>)iterator.getParent()).getIteratorIndex();
			
			SAMRecordList ll = mateBuffer.get(BamUtils.getPairId(rec));
			SAMRecord mate = findAndRemoveMate(rec, file,mateBuffer,ll);
			if (mate==null) {
//				if (rec.getMateAlignmentStart()-1>=pos) //otherwise there is no chance to find it later!
				if (fillMateBuffer)
					mateBuffer.put(BamUtils.getPairId(rec), new SAMRecordList(rec,ll,file));
//				System.out.println("to Matebuffer");
				return;
			}
//			if (ll==null) mateBuffer.remove(BamUtils.getPairId(rec));
			
			if (isAnyStrandSpecific()){
				SAMRecord first = rec.getFirstOfPairFlag()?rec:mate;
				SAMRecord second = rec.getFirstOfPairFlag()?mate:rec;
				if (!BamUtils.isValidStrand(ref.getStrand(), strandness[file], first,second))
					return;
				rec = first;
				mate = second;
			}
				
			
			if (!ignoreVariations && bamChecker!=null) {
				FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,mate,cumNumCond, pairedEndHandling==PairedEndHandling.JoinMates,false,ignoreVariations,keepReadNames, barcode);
				
				if (fac.isConsistent()) {
					fac.add(rec, mate, file);
					bamChecker.check(rec,mate,new ImmutableReferenceGenomicRegion<>(ref, fac,fac.create()));
				} else
					bamChecker.addInconsistent(rec,mate);
			}
			
			FactoryGenomicRegion fac = BamUtils.getFactoryGenomicRegion(rec,mate,cumNumCond, pairedEndHandling==PairedEndHandling.JoinMates,false,ignoreVariations,keepReadNames, barcode);
//			System.out.println("use it as "+fac.toRegionString());
			if (!fac.isConsistent()) {
//				System.out.println(rec.getSAMString());
//				System.out.println(mate.getSAMString());
//				IntArrayList coords1 = BamUtils.getGenomicRegionCoordinates(rec);
//				IntArrayList coords2 = BamUtils.getGenomicRegionCoordinates(mate);
//				GenomicRegion re1 = new ArrayGenomicRegion(coords1);
//				GenomicRegion re2 = new ArrayGenomicRegion(coords2);
//				System.out.println(re1+" "+re2+" "+re1.union(re2));
//				System.out.println();
				stats[4]++;
				return;
			}

			if (!map.containsKey(fac)) map.put(fac,fac);
			else fac = map.get(fac);
			
			stats[5]++;
			fac.add(rec, mate, file);
		}

		private SAMRecord findAndRemoveMate(SAMRecord rec, int file, 
				HashMap<String, SAMRecordList> mateBuffer, SAMRecordList ll) {
			if (ll==null) return null;
			if (BamUtils.checkMates(ll.r,rec) && file==ll.file) {
				if (ll.next==null) mateBuffer.remove(BamUtils.getPairId(rec));
				else mateBuffer.put(BamUtils.getPairId(rec),ll.next);
				return ll.r;
			}
			SAMRecordList last = ll;
			for (ll=ll.next; ll!=null; last = ll, ll=ll.next) {
				if (BamUtils.checkMates(ll.r,rec) && file==ll.file) {
					last.next = ll.next;
					return ll.r;
				}
			}
			return null;
		}

		private SAMRecord findAndRemoveMate(SAMRecord a, LinkedList<SAMRecord> ll) {
			Iterator<SAMRecord> it = ll.iterator();
			while (it.hasNext()) {
				SAMRecord c = it.next();
				if (a.getMateAlignmentStart()==c.getAlignmentStart()) {
					it.remove();
					return c;
				}
			}
			return null;
		}

		@Override
		public Spliterator<GenomicRegion> trySplit() {
//			if (!iterator().hasNext()) return null;
//			if (end-pos<minSpanForSplit) return null;
//
//			int center = (pos+end)/2;
//			if (center==pos) return null;
//
//			GenomicRegionSpliterator re = new GenomicRegionSpliterator(ref,center,end,false);
//			re.regionsReturned = regionsReturned;
//			re.basesCovered = basesCovered;
//			end = center;
//			return re;
			return null;
		}

		@Override
		public long estimateSize() {
			if (regionsReturned<START_ESTIMATE) return Long.MAX_VALUE;
			return (long) (Math.abs(end-pos)*(regionsReturned/(double)basesCovered));
		}

		@Override
		public int characteristics() {
			return DISTINCT|NONNULL|SORTED|ORDERED|IMMUTABLE;
		}

	}

	public int getNumConditions() {
		return cumNumCond[cumNumCond.length-1];
	}


	
	@Override
	public Spliterator<MutableReferenceGenomicRegion<AlignedReadsData>> iterateMutableReferenceGenomicRegions(
			ReferenceSequence reference) {

		Supplier<Function<GenomicRegion,MutableReferenceGenomicRegion<AlignedReadsData>>> supp = 
				()->{
					MutableReferenceGenomicRegion<AlignedReadsData> mut = new MutableReferenceGenomicRegion<AlignedReadsData>();
					return (e)->mut.set(reference, e, getAlignedReadsData(reference, e));
				};

				return new MappedSpliterator<GenomicRegion, MutableReferenceGenomicRegion<AlignedReadsData>>(iterateGenomicRegions(reference),supp);
	}

	@Override
	public Spliterator<MutableReferenceGenomicRegion<AlignedReadsData>> iterateIntersectingMutableReferenceGenomicRegions(
			ReferenceSequence reference, GenomicRegion region) {
		if (region.isEmpty()) return Spliterators.emptySpliterator();
		Supplier<Function<GenomicRegion,MutableReferenceGenomicRegion<AlignedReadsData>>> supp = 
				()->{
					MutableReferenceGenomicRegion<AlignedReadsData> mut = new MutableReferenceGenomicRegion<AlignedReadsData>();
					return (e)->mut.set(!isAnyStrandSpecific()?reference.toStrandIndependent():reference, e, getAlignedReadsData(reference, e));
				};

				return new MappedSpliterator<GenomicRegion, MutableReferenceGenomicRegion<AlignedReadsData>>(
						new FilteredSpliterator<GenomicRegion>(iterateGenomicRegions(reference, region.getStart(), region.getEnd()),r->region.intersects(r)),
						supp);
	}

	@Override
	public boolean add(ReferenceSequence reference, GenomicRegion region,
			AlignedReadsData data) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean remove(ReferenceSequence reference, GenomicRegion region) {
		throw new UnsupportedOperationException();
	}

	@Override
	public AlignedReadsData getData(ReferenceSequence reference,
			GenomicRegion region) {
		return getAlignedReadsData(reference, region);
	}
	
	@Override
	public long size() {
		return -1;
	}


	@Override
	public long size(ReferenceSequence reference) {
		MutableLong re = new MutableLong();
		iterateGenomicRegions(reference).forEachRemaining(g->re.N++);
		return re.N;
	}

	@Override
	public void clear() {
		throw new UnsupportedOperationException();
	}


	@Override
	public String toString() {
		return getName();
	}
	
	
}
