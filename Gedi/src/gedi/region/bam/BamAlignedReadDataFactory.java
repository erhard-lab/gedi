package gedi.region.bam;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.function.BiFunction;

import gedi.bam.tools.BamUtils;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.DefaultAlignedReadsData;
import gedi.core.region.GenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.sequence.DnaSequence;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;

public class BamAlignedReadDataFactory extends AlignedReadsDataFactory {

	private CharSequence genomicSequence;
	private HashMap<String,Integer> map = new HashMap<String, Integer>();
	private ArrayList<IntArrayList> distinctToIds = new ArrayList<IntArrayList>();
	private int[] cumNumCond;
	private GenomicRegion region;
	private boolean ignoreVariation;
	private boolean needReadNames;
	
	private int minIntronLength = -1;
	private boolean join;
	private boolean overlapping;
	
	BiFunction<SAMRecord,SAMRecord,DnaSequence> barcodeFun;
	
	public BamAlignedReadDataFactory(GenomicRegion region, int[] cumNumCond, boolean ignoreVariation, boolean needReadNames, boolean join, BiFunction<SAMRecord, SAMRecord, DnaSequence> barcode) {
		super(cumNumCond[cumNumCond.length-1]);
		this.cumNumCond = cumNumCond;
		this.region = region;
		this.ignoreVariation = ignoreVariation;
		this.needReadNames = needReadNames;
		this.join = join;
		this.barcodeFun = barcode;
	}
	
	public void setUseBlocks(int minIntronLength) {
		if (minIntronLength>=0 && !ignoreVariation) throw new RuntimeException("Cannot use blocks when recording variations!");
		this.minIntronLength = minIntronLength;
	}
	
	public void setBarcodeFun(BiFunction<SAMRecord,SAMRecord, DnaSequence> barcodeFun) {
		this.barcodeFun = barcodeFun;
	}
	
	@Override
	public AlignedReadsDataFactory start() {
		genomicSequence = null;
		map.clear();
		return super.start();
	}
	
	public boolean isOverlapping() {
		return overlapping;
	}
	
	public BamAlignedReadDataFactory setGenomicSequence(CharSequence genomicSequence) {
		this.genomicSequence = genomicSequence;
		return this;
	}
	
	public IntArrayList getReadIds(int distinct) {
		IntArrayList re = distinctToIds.get(distinct);
		re.sort();
		re.unique();
		return re;
	}
	
	public void addRecord(SAMRecord record, int file) {
		if (getCoveredGenomicLength(record)!=region.getTotalLength())
			throw new RuntimeException("Record and region do not match!");

		String key = getKey(record);
		Integer s = map.get(key);
		if (s==null) {
			// start a new distinct
			map.put(key, s = map.size());
			newDistinctSequence();
			
			Integer m = record.getIntegerAttribute("NH");
			if (m!=null) 
				setMultiplicity(s, m);
			else
				setMultiplicity(s, 0);
			
			if (!ignoreVariation) {
				boolean second = record.getReadPairedFlag() && record.getSecondOfPairFlag();
				addVariations(record,record.getReadNegativeStrandFlag()!=second,record.getReadNegativeStrandFlag(),0,getCurrentVariationBuffer());
			}
			else
				setMultiplicity(s, 0);
			
			if (needReadNames)
				distinctToIds.add(new IntArrayList());
		}
		
		if (needReadNames) {
//			if (record.getReadPairedFlag() && record.getFirstOfPairFlag())
//				distinctToIds.get(s).add(record.getReadName()+BamUtils.FIRST_PAIR_SUFFIX);
//			else if (record.getReadPairedFlag() && record.getSecondOfPairFlag())
//				distinctToIds.get(s).add(record.getReadName()+BamUtils.SECOND_PAIR_SUFFIX);
//			else 
			if (!StringUtils.isInt(StringUtils.splitField(record.getReadName(),'#',0)))
				throw new RuntimeException("Can only keep integer read names!");
			
			distinctToIds.get(s).add(Integer.parseInt(StringUtils.splitField(record.getReadName(),'#',0)));
			
		}

		// handle counts
		int start = file==0?0:cumNumCond[file-1];
		int end = cumNumCond[file];
		
		String c = record.getAttribute("XR") instanceof String?record.getStringAttribute("XR"):null;
		if (c==null && end-start!=1) throw new RuntimeException("XR counts do not match BAM header descriptor");
		if (c==null)
			incrementCountSingle(s, start, barcodeFun!=null?barcodeFun.apply(record, null):null);
		else {
			String[] f = StringUtils.split(c,',');
			if (end-start!=f.length) throw new RuntimeException("XR counts do not match BAM header descriptor");
			for (int i=0; i<f.length; i++)
				incrementCount(s, start+i, Integer.parseInt(f[i]));
		}
		
		
	}

	@Override
	public DefaultAlignedReadsData create() {
		if (needReadNames) {
			for (int d=0; d<distinctToIds.size(); d++)
				setId(d, getReadIds(d).getInt(0));
		}
		return super.create();
	}
	
	private String getKey(SAMRecord record) {
		if (ignoreVariation) return "SAM";
		
		String xv = record.getStringAttribute("XV");
		String secondary = xv==null?record.getReadString():xv;
		
		return record.getCigarString()+secondary;
	}

	/**
	 * For mate pairs; The region of this must match the two records: intron consistent, contained, start with record1 and end with record2
	 * 
	 * If mappings overlap, variations are resolved optimistically (i.e. mismatch only in one mate pair -> to variation); if there are are 
	 * inconsistent variations, the variation of the first mate is taken.
	 * @param record1
	 * @param record2
	 * @param file
	 */
	public void addRecord(SAMRecord record1, SAMRecord record2, int file) {
//		if (record1.getAlignmentStart()>record2.getAlignmentStart() ||
//				(record1.getAlignmentStart()==record2.getAlignmentStart() && record1.getAlignmentEnd()>record2.getAlignmentEnd())
//				) {
//			SAMRecord d = record2;
//			record2 = record1;
//			record1 = d;
//		}
		// this would break the calculation of the mismatch positions; this was here for the key (now handled below!)
		
		GenomicRegion reg1 = minIntronLength>=0?BamUtils.getArrayGenomicRegionByBlocks(record1,minIntronLength):BamUtils.getArrayGenomicRegion(record1);
		GenomicRegion reg2 = minIntronLength>=0?BamUtils.getArrayGenomicRegionByBlocks(record2,minIntronLength):BamUtils.getArrayGenomicRegion(record2);
		
		if (!region.isIntronConsistent(reg1) || !region.isIntronConsistent(reg2) || !reg1.isIntronConsistent(reg2) 
				|| !region.contains(reg1) || !region.contains(reg2) ) {
//				|| region.induce(reg1.getStart())!=0 || region.induce(reg2.getStop())!=region.getTotalLength()-1)
			// can be for softclipped reads
			
			// Stupid STAR: this happens!!! 4776729-4776803 and 4776750-4776801|4777524-4777549 are pairs
			throw new RuntimeException("Record and regions do not match: \n"
					+ "Region: "+region+"\nRecord1: "+reg1+"\t"+record1.getSAMString()
					+"Record2: "+reg2+"\t"+record2.getSAMString());
		}
		
		String key = getKey(record1)+getKey(record2);
		if (record1.getAlignmentStart()>record2.getAlignmentStart() ||
				(record1.getAlignmentStart()==record2.getAlignmentStart() && record1.getAlignmentEnd()>record2.getAlignmentEnd())
				) {
			key = getKey(record2)+getKey(record1);
		}
		
		Integer s = map.get(key);
		if (s==null) {
			// start a new distinct
			map.put(key, s = map.size());
			newDistinctSequence();
			if (needReadNames)
				distinctToIds.add(new IntArrayList());
			
			Integer m = record1.getIntegerAttribute("NH");
			if (m!=null) 
				setMultiplicity(s, m);
			else
				setMultiplicity(s, 0);
			
			if (!join)
				setGeometry(s, reg1.subtract(reg2).getTotalLength(), reg1.intersect(reg2).getTotalLength(), reg2.subtract(reg1).getTotalLength());
			
			if (reg1.intersects(reg2)) {
				overlapping = true;
//				GenomicRegion off = region.induce(reg1.subtract(reg2));
				int off1, off2;
				if (!record1.getReadNegativeStrandFlag()) {
					off1 = region.induce(reg1.getStart());
					off2 = region.induce(reg2.getStart());
				} else {
					off1 = region.getTotalLength()-1-region.induce(reg1.getStop());
					off2 = region.getTotalLength()-1-region.induce(reg2.getStop());
				}

				if (!ignoreVariation) {
//					addVariations(record1,record1.getReadNegativeStrandFlag(),record1.getReadNegativeStrandFlag(),0,getCurrentVariationBuffer());
//					addVariations(record2,record1.getReadNegativeStrandFlag(),record2.getReadNegativeStrandFlag(),off.getTotalLength(),getCurrentVariationBuffer());
					addVariations(record1,record1.getReadNegativeStrandFlag(),record1.getReadNegativeStrandFlag(),off1,getCurrentVariationBuffer());
					addVariations(record2,record1.getReadNegativeStrandFlag(),record2.getReadNegativeStrandFlag(),off2,getCurrentVariationBuffer());
					
					// to merge variations in the overlap;
					/*
					ArrayList<VarIndel> buffer1 = new ArrayList<VarIndel>();
					ArrayList<VarIndel> buffer2 = new ArrayList<VarIndel>();
					addVariations(record1,record1.getReadNegativeStrandFlag(),record1.getReadNegativeStrandFlag(),0,buffer1);
					addVariations(record2,record1.getReadNegativeStrandFlag(),record2.getReadNegativeStrandFlag(),off.getTotalLength(),buffer2);
					buffer1.sort(FunctorUtils.naturalComparator());
					buffer2.sort(FunctorUtils.naturalComparator());
					
					for (int i1=0, i2=0; i1<buffer1.size() || i2<buffer2.size(); ) {
						int p1 = i1<buffer1.size()?buffer1.get(i1).getPosition():Integer.MAX_VALUE;
						int p2 = i2<buffer2.size()?buffer2.get(i2).getPosition():Integer.MAX_VALUE;
						
						if (p1==p2) {
							getCurrentVariationBuffer().add(buffer1.get(i1));
							i1++; i2++;
						}
						else if (p1<p2) {
							if (!overlap.contains(p1))
								getCurrentVariationBuffer().add(buffer1.get(i1));
							i1++;
						}
						else {
							if (!overlap.contains(p2))
								getCurrentVariationBuffer().add(buffer2.get(i2));
							i2++;
						}
					}
					*/
				}
				else
					setMultiplicity(s, 0);
			} else {
				if (!ignoreVariation) {
					addVariations(record1,record1.getReadNegativeStrandFlag(),record1.getReadNegativeStrandFlag(),0,getCurrentVariationBuffer());
					addVariations(record2,record1.getReadNegativeStrandFlag(),record2.getReadNegativeStrandFlag(),reg1.getTotalLength(),getCurrentVariationBuffer());
				}
				else
					setMultiplicity(s, 0);
			}
			
		}

		if (!BamUtils.getPairId(record1).equals(BamUtils.getPairId(record2)))
				throw new RuntimeException("Read names do not match for mate pairs!");
		if (needReadNames) {
			if (!StringUtils.isInt(StringUtils.splitField(record1.getReadName(),'#',0)))
				throw new RuntimeException("Can only keep integer read names!");
			distinctToIds.get(s).add(Integer.parseInt(StringUtils.splitField(record1.getReadName(),'#',0)));
		}
		
		// handle counts
		int start = file==0?0:cumNumCond[file-1];
		int end = cumNumCond[file];
		
		String c = record1.getAttribute("XR") instanceof String?record1.getStringAttribute("XR"):null;
		if (c==null && end-start!=1) throw new RuntimeException("XR counts to not match BAM header descriptor");
		if (c==null)
			incrementCountSingle(s, start,barcodeFun!=null?barcodeFun.apply(record1, record2):null);
		else {
			String[] f = StringUtils.split(c,',');
			if (end-start!=f.length) throw new RuntimeException("XR counts to not match BAM header descriptor");
			for (int i=0; i<f.length; i++)
				incrementCount(s, start+i, Integer.parseInt(f[i]));
		}
		
	}
	
	
	/**
	 * Offset is for paired end reads to specify the position, where thes second records starts within the genomic region.
	 * @param record
	 * @param offset
	 */
	private void addVariations(SAMRecord record, boolean invertPos, boolean complement, int offset, Collection<VarIndel> to) {
		
		String xv = record.getStringAttribute("XV");
		if (xv!=null) {
			for (String v : StringUtils.split(xv, ','))
				to.add(createVarIndel(v));
			return;
		}
		
		if (record.getStringAttribute("MD")==null || record.getReadString().equals("*")) 
			throw new RuntimeException("Need read sequences and MD tag to read variations for reads!");
		
		int coveredGenomic = getCoveredGenomicLength(record);
		
		int pos;
		CharSequence vread;
		
		int ls = 0;
		int lr = 0;
		for (CigarElement e : record.getCigar().getCigarElements()) {
			switch (e.getOperator()){
			case I: 
				pos = (invertPos?coveredGenomic-lr:lr);
				vread = complement?SequenceUtils.getDnaComplement(getReadSequence(record,ls,ls+e.getLength())):getReadSequence(record,ls,ls+e.getLength());;
				if (invertPos)
					vread = StringUtils.reverse(vread);
				ls+=e.getLength(); 
				to.add(createInsertion(pos+offset, vread, record.getReadPairedFlag() && record.getSecondOfPairFlag()));
				break;
			case D:
				pos = (invertPos?coveredGenomic-lr-e.getLength():lr);
				vread = complement?SequenceUtils.getDnaReverseComplement(getReferenceSequence(record,lr,lr+e.getLength())):getReferenceSequence(record,lr,lr+e.getLength());;
				to.add(createDeletion(pos+offset, vread, record.getReadPairedFlag() && record.getSecondOfPairFlag()));
				lr+=e.getLength(); 
				break;
			case M: 
				CharSequence read = getReadSequence(record,ls,ls+e.getLength());
				CharSequence ref = getReferenceSequence(record,lr,lr+e.getLength());
				for (int i=0; i<read.length(); i++)
					if (read.charAt(i)!=ref.charAt(i)) {
						pos = (invertPos?coveredGenomic-1-lr-i:lr+i);
						char g = complement?SequenceUtils.getDnaComplement(ref.charAt(i)):ref.charAt(i);
						char r = complement?SequenceUtils.getDnaComplement(read.charAt(i)):read.charAt(i);
						to.add(createMismatch(pos+offset, g, r, record.getReadPairedFlag() && record.getSecondOfPairFlag()));
					}
				ls+=e.getLength();
				lr+=e.getLength();
				break;
			case S:
				pos = (invertPos?coveredGenomic-lr:lr);
				vread = complement?SequenceUtils.getDnaComplement(getReadSequence(record,ls,ls+e.getLength())):getReadSequence(record,ls,ls+e.getLength());;
				if (invertPos)
					vread = StringUtils.reverse(vread);
				to.add(createSoftclip(pos==0, vread, record.getReadPairedFlag() && record.getSecondOfPairFlag()));
				
				ls+=e.getLength();
				
				break;
			case N: break;
			default: throw new IllegalArgumentException("Cigar operator "+e.getOperator()+" unknown!");
			}
		}
		
	}


	private int getCoveredGenomicLength(SAMRecord record) {
		int re = 0;
		for (CigarElement e : record.getCigar().getCigarElements()) {
			switch (e.getOperator()){
			case D: 
			case M: 
				re+=e.getLength();
				break;
			case N:
			case S:
			case I: 
				break;
			default: throw new IllegalArgumentException("Cigar operator "+e.getOperator()+" unknown!");
			}
		}
		return re;
	}

	private CharSequence getReferenceSequence(SAMRecord r, int s, int e) {
//		String md = r.getStringAttribute("MD");
//		if (md==null || r.getReadString().equals("*")) {
//			
//			
//			if (genomicSequence==null) {
//				if (md!=null) {
//					MismatchString mm = new MismatchString(r.getStringAttribute("MD"));
//					String seq = StringUtils.repeatSequence('N',mm.getGenomicLength()).toString();
//					seq = mm.reconstitute(seq);
//					return seq.substring(s, e);
//				}
//				return StringUtils.repeatSequence('N', e-s);
//			}
//			return genomicSequence.subSequence(s, e);
//		}
		return BamUtils.restoreSequence(r,false).substring(s, e);
	}
	private CharSequence getReadSequence(SAMRecord r, int s, int e) {
//		if (r.getReadString().length()==0 || r.getReadString().equals("*")) return StringUtils.repeatSequence('N', e-s);
		return r.getReadString().substring(s,e);
	}
	
	
	
}
