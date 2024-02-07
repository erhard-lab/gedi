package gedi.slam;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.function.Predicate;

import gedi.core.data.annotation.NameProvider;
import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Strandness;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.mutable.MutableDouble;

public class SlamCollector {

	private final static String[] nil = new String[0];
	private Genomic genomic;
	private GenomicRegionStorage<AlignedReadsData> reads;
	private GenomicRegionStorage<?> masked;
	private GenomicRegionStorage<NameProvider> locations;
	
	private HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<Transcript>>> gene2Trans;
	private int cond;
	
	private ArrayList<MismatchCounter> counter = new ArrayList<>();
	private MismatchCounter modelCounter;
	private MismatchCounter objectCounter;
	private MismatchCounter intronCounter;
	
	private Strandness strandness;
	private AtomicInteger readLength1 = new AtomicInteger(0);
	private AtomicInteger readLength2 = new AtomicInteger(0);
	
	private int trim5p=25;
	private int trim3p=0;
	
	private ReadCountMode mode,overlap;
	
	private int verbose = 0;
	private BiConsumer<ImmutableReferenceGenomicRegion<String>, int[]> numiReport;
	private boolean[] no4sU;
	
	private boolean lenientOverlap;
	
	private boolean useAllReadsForModel = false;
	private boolean highmem = false;
	
	public SlamCollector(Genomic genomic, Predicate<String> keepGene, GenomicRegionStorage<AlignedReadsData> reads, GenomicRegionStorage<?> masked, GenomicRegionStorage<NameProvider> locations, Strandness strandness, int trim5p, int trim3p, ReadCountMode mode, ReadCountMode overlap, boolean[] no4sU, boolean countIntrons, boolean lenientOverlap, boolean modelall, boolean highmem) {
		this.genomic = genomic;
		this.reads = reads;
		this.masked = masked;
		this.locations = locations;
		this.strandness = strandness;
		this.trim5p = trim5p;
		this.trim3p = trim3p;
		this.mode = mode;
		this.overlap = overlap;
		this.no4sU = no4sU;
		this.lenientOverlap = lenientOverlap;
		this.useAllReadsForModel = modelall;
		this.highmem = highmem;
		
		this.cond = reads.getMetaDataConditions().length;
		gene2Trans = genomic.getTranscripts().ei().filter(r->keepGene.test(r.getData().getGeneId())).indexMulti(t->t.getData().getGeneId(), t->t);
		
		counter.add(objectCounter = new MismatchCounter("Exonic", cond, r->r.genes.length>0, r->r.genes));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) 
			counter.add(new MismatchCounter("ExonicSense", cond, r->r.genes.length>0 && !r.isOppositeStrand(), r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) 
			counter.add(new MismatchCounter("ExonicAntisense", cond, r->r.genes.length>0 && r.isOppositeStrand(), r->nil));
		if (locations!=null) {
			counter.add(objectCounter = new MismatchCounter("Locations", cond, r->r.locations.length>0, r->r.locations));
			if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) 
				counter.add(new MismatchCounter("LocationsSense", cond, r->r.locations.length>0 && !r.isOppositeStrand(), r->nil));
			if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) 
				counter.add(new MismatchCounter("LocationsAntisense", cond, r->r.locations.length>0 && r.isOppositeStrand(), r->nil));
		}
		counter.add(intronCounter = new MismatchCounter("Intronic", cond, r->
								r.locations.length==0 
								&& r.genes.length==0 
								&& r.read.getRegion().getNumParts()==1
								&& EI.wrap(r.tr).filter(t->t.getRegion().invert().intersects(r.read.getRegion())).count()>1,
							r->EI.wrap(r.overlapgenes).map(n->n+"_intronic").toArray(String.class)));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) 
			counter.add(new MismatchCounter("IntronicSense", cond, r->
								r.locations.length==0 
								&& r.genes.length==0 
								&& r.read.getRegion().getNumParts()==1
								&& EI.wrap(r.tr).filter(t->t.getRegion().invert().intersects(r.read.getRegion())).count()>1
								&& !r.isOppositeStrand(),
							r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) 
			counter.add(new MismatchCounter("IntronicAntisense", cond, r->
								r.locations.length==0 
								&& r.genes.length==0 
								&& r.read.getRegion().getNumParts()==1
								&& EI.wrap(r.tr).filter(t->t.getRegion().invert().intersects(r.read.getRegion())).count()>1
								&& r.isOppositeStrand(),
							r->nil));
		
		if (!countIntrons)
			intronCounter = null;
		
		modelCounter = objectCounter;
		if (modelall) {
			counter.add(modelCounter = new MismatchCounter("All", cond, r->true,r->nil));
		}
		
	}
	
	public void countNumis(BiConsumer<ImmutableReferenceGenomicRegion<String>, int[]> report) {
		numiReport = report;
	}
	
	public Collection<String> getCounterTypes() {
		return EI.wrap(counter).map(m->m.name).list();
	}

	
	public void ercc() {
		counter.add(new MismatchCounter("ERCC", cond, r->r.read.getReference().getName().startsWith("ERCC"),r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) 
			counter.add(new MismatchCounter("ERCCSense", cond, r->r.read.getReference().getName().startsWith("ERCC") && !r.isOppositeStrand(),r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) 
			counter.add(new MismatchCounter("ERCCAntisense", cond, r->r.read.getReference().getName().startsWith("ERCC") && r.isOppositeStrand(),r->nil));
	}

	
	public void addReferenceCounter(String label, String name, boolean invert) {
		counter.add(new MismatchCounter(label, cond, r->r.genes.length>0 && (r.read.getReference().getName().equals(name)^invert),r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) 
			counter.add(new MismatchCounter(label+"Sense", cond, r->r.genes.length>0 && (r.read.getReference().getName().equals(name)^invert) && !r.isOppositeStrand(),r->nil));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) 
			counter.add(new MismatchCounter(label+"Antisense", cond, r->r.genes.length>0 && (r.read.getReference().getName().equals(name)^invert) && r.isOppositeStrand(),r->nil));
	}

	public void setVerbose(int verbose) {
		this.verbose = verbose;
	}
	
	public int getTotalReadLength() {
		return readLength1.get()+readLength2.get();
	}
	
	public ExtendedIterator<GeneData> collect(String gene){
		ReferenceGenomicRegion<Transcript> rr = (ReferenceGenomicRegion<Transcript>) genomic.getNameIndex().get(gene);
		if (rr==null)
			rr = (ReferenceGenomicRegion<Transcript>) genomic.getNameIndex().getUniqueWithPrefix(gene);
		if (rr==null)
			rr = ImmutableReferenceGenomicRegion.parse(gene);
		ReferenceGenomicRegion<Transcript> rru = rr; 
//		return collect(genomic.getGenes().ei(rr).filter(x->x.compareTo(rru)==0).getUniqueResult(true, true));		
		return collect(new ImmutableReferenceGenomicRegion<>(rr.getReference(), rr.getRegion(), rr.getData().getGeneId()).toImmutable());		
	}
	public ExtendedIterator<GeneData> collect(ImmutableReferenceGenomicRegion<String> gene){
		if (locations!=null) { 
			GenomicRegion union = locations.ei(gene).reduce(gene.getRegion(),(l,u)->u.union(l.getRegion()));
			gene = new ImmutableReferenceGenomicRegion<>(gene.getReference(), union, gene.getData());
		}
		ArrayList<ImmutableReferenceGenomicRegion<Transcript>> tr = gene2Trans.get(gene.getData());
		
		
		// create local copies for thread safety
		MismatchCounter[] counter = new MismatchCounter[this.counter.size()];
		for (int i=0; i<counter.length; i++)
			counter[i] = this.counter.get(i).createForThread();
		int readLen1 = readLength1.get();
		int readLen2 = readLength2.get();
		
		NumiAlgorithm numi = numiReport==null?null:new NumiAlgorithm(gene);
		
		
		// buffer for all counters of the current read
		ArrayList<MismatchCounter> useCounters = new ArrayList<>();
		
		// Map to store mismatches for the objects to count (transcripts or locations)
		HashMap<String,MismatchCounter> objectCounter = new HashMap<>();
		
		char[] seq = genomic.getSequence(gene).toString().toUpperCase().toCharArray();
		IntPair pp = new IntPair(-1, -1);
		ReadInfo readinfo = new ReadInfo(gene,tr);
		
		MutableSingletonGenomicRegion checkerReg = new MutableSingletonGenomicRegion();
		MutableReferenceGenomicRegion<Void> checker = new MutableReferenceGenomicRegion<Void>()
						.setReference(gene.getReference().toStrandIndependent())
						.setRegion(checkerReg);
		
		char[] overlapfound = new char[2000]; // reads arent longer than this! read sense base
		
		ExtendedIterator<ImmutableReferenceGenomicRegion<AlignedReadsData>> it = EI.empty();
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Sense)) it = it.chain(reads.ei(gene));
		if (strandness.equals(Strandness.Unspecific) || strandness.equals(Strandness.Antisense)) it = it.chain(reads.ei(gene.toMutable().toOppositeStrand()));
		
		ReadCountMode mode = this.mode;
		
		if (highmem) it = EI.wrap(it.list());
		
		for (ImmutableReferenceGenomicRegion<AlignedReadsData> read : it.loop()) {
//			if (read.toLocationString().contains("70|28260"))
//					System.out.println();
			useCounters.clear();

			// set readinfo
			readinfo.read = read;
			// tr only contains transcripts of this gene, so this is either empty or contains this gene!
			readinfo.genes = EI.wrap(tr).filter(t->isConsistent(t,read)).map(t->t.getData().getGeneId()).unique(false).toArray(String.class);
			readinfo.overlapgenes = EI.wrap(tr).map(t->t.getData().getGeneId()).unique(false).toArray(String.class);
			
			if (overlap!=ReadCountMode.All) {
				if (overlap==ReadCountMode.CollapseAll) {
					mode = this.mode.transformCounts(c->1.0);
				}
				else {
					HashSet<String> gset = genomic.getTranscripts().ei(read).filter(t->isConsistent(t, read)).map(t->t.getData().getGeneId()).set();
					if (gset.size()>1) {
						if (overlap==ReadCountMode.Unique || overlap==ReadCountMode.CollapseUnique) continue;
						if (overlap==ReadCountMode.Divide || overlap==ReadCountMode.Weight) {
							mode = this.mode.transformCounts(c->c/gset.size());
						}
					} else if (overlap==ReadCountMode.CollapseUnique)
						mode = this.mode.transformCounts(c->1.0);
				}
			}
			
			if (locations!=null)
				readinfo.locations = locations.ei(gene).filter(l->l.getRegion().contains(read.getRegion())).map(r->r.getData().getName()).unique(false).toArray(String.class);
			
			for (MismatchCounter c : counter)
				if (c.predicate.test(readinfo))
					useCounters.add(c);

			if (!useCounters.isEmpty()) {
				AlignedReadsData rd = read.getData();
				char conversionBase = 'T';
				char[] buff = gene.getRegion().contains(read.getRegion())?SequenceUtils.extractSequence(gene.induce(read.getRegion()), seq):genomic.getSequenceSave(read.toMutable().toStrand(gene.getReference().getStrand())).toString().toUpperCase().toCharArray();
				
				if (readinfo.isOppositeStrand()) {
					buff = SequenceUtils.getDnaReverseComplement(new String(buff)).toCharArray();
					conversionBase = 'A';
				}
					
				int n = read.getRegion().getTotalLength();
				
				for (int d=0; d<rd .getDistinctSequences(); d++) {
					
					int conv = 0;
					int doubleconv = 0;
					
					int trl1 = read.getData().hasGeometry()?read.getData().getReadLength1(d):read.getRegion().getTotalLength();
					if (trl1>readLen1) 
						readLen1=this.readLength1.accumulateAndGet(trl1, Math::max);
					
					int trl2 = read.getData().hasGeometry()?read.getData().getReadLength2(d):read.getRegion().getTotalLength();
					if (trl2>readLen2) 
						readLen2=this.readLength2.accumulateAndGet(trl2, Math::max);
					
					// we may make a small error if the first read(s) do not have maximal length, but that is not too bad!
					
					
					for (int v=0; v<rd.getVariationCount(d); v++) {
						
						if (rd.isMismatch(d, v)) {
							
							if (!isBase(rd.getMismatchGenomic(d, v).charAt(0)) || !isBase(rd.getMismatchRead(d, v).charAt(0)))
								continue;
							
							checkerReg.p = read.map(rd.getMismatchPos(d, v));
							if (masked==null || masked.ei(checker).count()==0) {
								int pos = rd.getMismatchPos(d, v);
								int mpos = rd.mapToRead(d, rd.getMismatchPos(d, v), rd.isVariationFromSecondRead(d, v), readLen1+readLen2);
								if (!rd.isPositionInOverlap(d, pos)) {
									if (mpos<trim5p || (read.getData().hasGeometry() && mpos>=readLen1+readLen2-trim5p))
										continue;
									if ((readLen1-mpos<trim3p && mpos<readLen1) || (read.getData().hasGeometry() && mpos>=readLen1 && mpos<readLen1+trim3p))
										continue;
								}
								
								char checkGenomic = 'T';
								char checkRead = 'C';
								// this is necessary here, as the check does not look at buff, but in the ard
								if (!readinfo.isOppositeStrand()) {
									if (rd.isVariationFromSecondRead(d, v)) {
										checkGenomic = 'A';
										checkRead = 'G';
									}
								} else {
									if (!rd.isVariationFromSecondRead(d, v)) {
										checkGenomic = 'A';
										checkRead = 'G';
									}
								}
								boolean thisconv = false;
								boolean thisdoubleconv = false;
								if (rd.getMismatchGenomic(d, v).charAt(0)==checkGenomic  && rd.getMismatchRead(d, v).charAt(0)==checkRead && buff[pos]!='N') {
									if (rd.isPositionInOverlap(d, pos)) {
										if (overlapfound[pos]=='C') { // double hit: set save
											doubleconv++;
											thisdoubleconv = true;
											if (numi!=null && readinfo.genes.length>0 && rd.getMultiplicity(d)<=1) numi.addConversion(read.map(rd.getMismatchPos(d, v)));
										}
									} else {
										conv++;
										thisconv = true;
									}
								}
								
								if (rd.isPositionInOverlap(d, pos)) {
									char genomicBase  = rd.getMismatchGenomic(d, v).charAt(0);
									char readBase = rd.getMismatchRead(d, v).charAt(0);
									if (readinfo.isOppositeStrand() ^ rd.isVariationFromSecondRead(d, v)) {
										readBase = SequenceUtils.getDnaComplement(readBase);
										genomicBase = SequenceUtils.getDnaComplement(genomicBase);
									}
									if (overlapfound[pos]==readBase)
										for (MismatchCounter c : useCounters)
											c.countDoubleHit(SequenceUtils.inv_nucleotides[genomicBase],SequenceUtils.inv_nucleotides[readBase],rd,d, mode);
									overlapfound[pos]=readBase;
								}
								
								
								if ((!rd.isVariationFromSecondRead(d, v) && Character.toUpperCase(buff[rd.getMismatchPos(d, v)])!=Character.toUpperCase(rd.getMismatchGenomic(d, v).charAt(0)) && Character.toUpperCase(buff[rd.getMismatchPos(d, v)])!='N')
										|| (rd.isVariationFromSecondRead(d, v) && Character.toUpperCase(buff[rd.getMismatchPos(d, v)])!=Character.toUpperCase(SequenceUtils.getDnaComplement(rd.getMismatchGenomic(d, v).charAt(0))) && Character.toUpperCase(buff[rd.getMismatchPos(d, v)])!='N')
										)
									throw new RuntimeException("Mismatch characters in mapping file are not consistent with genome: Genome="+new String(buff)+" "+rd.getVariation(d, v)+"\n"+read+"\n@ "+gene);
								
								int indg = SequenceUtils.inv_nucleotides[rd.getMismatchGenomic(d, v).charAt(0)];
								int indr = SequenceUtils.inv_nucleotides[rd.getMismatchRead(d, v).charAt(0)];
								pos = rd.mapToRead(d, rd.getMismatchPos(d, v), rd.isVariationFromSecondRead(d, v), readLen1+readLen2);
//								if (indg==2 && indr==0 && pos==151 && !readinfo.isOppositeStrand() && rd.isPositionInOverlap(d, rd.getMismatchPos(d, v)))
//									System.out.println(read+" "+rd.mapToRead(d, rd.getMismatchPos(d, v), rd.isVariationFromSecondRead(d, v), readLen1));
								if (indg!=indr && Math.min(indg, indr)>=0 && Math.max(indg, indr)<4) {
									if ((verbose>0 && (thisconv || thisdoubleconv)) || verbose>1) 
										System.out.println(read.getReference()+":"+read.map(rd.getMismatchPos(d, v))+"\t"+rd.getMismatchGenomic(d, v).charAt(0)+"->"+rd.getMismatchRead(d, v).charAt(0)+" "+(rd.isPositionInFirstRead(d, rd.getMismatchPos(d, v))?1:0)+(rd.isPositionInSecondRead(d, rd.getMismatchPos(d, v))?1:0)+"\t"+StringUtils.toString(readinfo.genes)+"\t"+(thisconv?1:0)+"\t"+(thisdoubleconv?1:0)+"\t"+rd.getCountsForDistinct(null,d, mode).toArrayString("\t", false));
									for (MismatchCounter c : useCounters)
										c.countMismatch(indg,indr,rd,d, pos,rd.isPositionInFirstRead(d, rd.getMismatchPos(d, v)),rd.isPositionInSecondRead(d, rd.getMismatchPos(d, v)),readinfo.isOppositeStrand(), mode, no4sU);
								}
							} 
//							else
//								if (rd.getMismatchGenomic(d, v).charAt(0)!='T')
//									throw new RuntimeException();
							// cannot do this anymore, all SNPs are identified (but this DID work!)
						}
					}
					
					if (numi!=null) {
						numi.finishRead(read.map(new ArrayGenomicRegion(rd.getGeometryBeforeOverlap(d),rd.getGeometryBeforeOverlap(d)+rd.getGeometryOverlap(d))),rd,d);
					}
					
					// clear overlap
					for (int v=0; v<rd.getVariationCount(d); v++) 
						if (rd.isMismatch(d, v)) 
							overlapfound[rd.getMismatchPos(d, v)]=0;
					
					int ttotal = 0;
					int ttotaldouble = 0;
					for (int i=0; i<n; i++) {
						
						int mpos = rd.mapToRead(d, i, false, readLen1+readLen2);
						if (!rd.isPositionInOverlap(d, i)) {
							if (mpos<trim5p || (rd.hasGeometry() && mpos>=readLen1+readLen2-trim5p))
								continue;
							if ((readLen1-mpos<trim3p && mpos<readLen1) || (read.getData().hasGeometry() && mpos>=readLen1 && mpos<readLen1+trim3p))
								continue;
						}
						
						
						// check if masked
						checkerReg.p = read.map(i);
						if ((masked==null || masked.ei(checker).count()==0) && isBase(buff[i])) {
							if (rd.isPositionInFirstRead(d, i)) {
								if (buff[i]==conversionBase && !rd.isPositionInOverlap(d, i))
									ttotal++;
								int genomic = SequenceUtils.inv_nucleotides[buff[i]];
								if (genomic>=0 && genomic<4) {
									
									int pos = rd.mapToRead1(d, i);
									for (MismatchCounter c : useCounters) {
										c.countTotal(genomic,rd,d, pos,rd.isPositionInFirstRead(d, i),rd.isPositionInSecondRead(d, i),readinfo.isOppositeStrand(), mode, no4sU);
//										c.countTotal(genomic, rd, d);
									}
								}
							}
							if (rd.isPositionInSecondRead(d, i)) {
								if (buff[i]==conversionBase && !rd.isPositionInOverlap(d, i))
									ttotal++;
								int genomic = SequenceUtils.inv_nucleotides[SequenceUtils.getDnaComplement(buff[i])];
								if (genomic>=0 && genomic<4) {
									int pos = rd.mapToRead2(d, i,readLen1+readLen2);
									for (MismatchCounter c : useCounters) {
										c.countTotal(genomic,rd,d, pos,rd.isPositionInFirstRead(d, i),rd.isPositionInSecondRead(d, i),readinfo.isOppositeStrand(), mode, no4sU);
//										c.countTotal(genomic, rd, d);
									}
								}
							}
							
							if (rd.isPositionInFirstRead(d, i) && rd.isPositionInSecondRead(d, i)) {
								if (buff[i]==conversionBase)
									ttotaldouble++;
								int ind = SequenceUtils.inv_nucleotides[buff[i]];
								if (ind>=0 && ind<4)
									for (MismatchCounter c : useCounters) 
										c.countDoubleTotal(ind, rd, d, mode);
							}
								
						}
					}
					if (conv>ttotal)
						throw new RuntimeException("Can never ever ever happen!"+read+" "+d);
					 
//					if (doubleconv>0) {
//						if (verbose && (conv>0 || doubleconv>0)) 
//							System.out.println(read+"\t"+conv+"\t"+doubleconv+"\t"+rd.getTotalCountsForConditions(null, mode).toArrayString("\t", false)+"\n");
						
//						if (doubleconv>0) {
//							pp.set(ttotal, conv);
//							for (MismatchCounter c : useCounters)
//								c.countOutsideWhenDouble(pp, rd, d, mode);
//							for (String object : this.objectCounter.objectsToCount.apply(readinfo))
//								objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countOutsideWhenDouble(pp, rd, d, mode);
//						}
						
						if (ttotaldouble>0) {
							if (verbose>2) 
								System.out.println(read+"\t"+StringUtils.toString(readinfo.genes)+"\t"+doubleconv+"\t"+ttotaldouble+"\tdouble");
							
							pp.set(ttotaldouble, doubleconv);
							for (MismatchCounter c : useCounters)
								c.countDoubleHitBinomial(pp, rd, d, mode);
							
							for (String object : this.objectCounter.objectsToCount.apply(readinfo))
								objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countDoubleHitBinomial(pp, rd, d, mode);
							
							if (this.intronCounter!=null && this.intronCounter.predicate.test(readinfo))
								for (String object : this.intronCounter.objectsToCount.apply(readinfo))
									objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countDoubleHitBinomial(pp, rd, d, mode);
						}
						
						pp.set(ttotaldouble+ttotal, doubleconv+conv);
						for (String object : this.objectCounter.objectsToCount.apply(readinfo))
							objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countBothBinomial(pp, rd, d, mode);
						
						if (this.intronCounter!=null && this.intronCounter.predicate.test(readinfo))
							for (String object : this.intronCounter.objectsToCount.apply(readinfo))
								objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countBothBinomial(pp, rd, d, mode);
//					}
//					else if (ttotal>0) {
//						if (verbose && (conv>0 || doubleconv>0)) 
//							System.out.println(read+"\t"+conv+"\t"+doubleconv+"\t"+rd.getTotalCountsForConditions(null, mode).toArrayString("\t", false)+"\n");
						
						if (verbose>2) 
							System.out.println(read+"\t"+StringUtils.toString(readinfo.genes)+"\t"+conv+"\t"+ttotal+"\tsingle");
						
						pp.set(ttotal, conv);
						for (MismatchCounter c : useCounters)
							c.countBinomial(pp, rd, d, mode);
						
						for (String object : this.objectCounter.objectsToCount.apply(readinfo))
							objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countBinomial(pp, rd, d, mode);
						
						if (this.intronCounter!=null && this.intronCounter.predicate.test(readinfo))
							for (String object : this.intronCounter.objectsToCount.apply(readinfo))
								objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).countBinomial(pp, rd, d, mode);
//					}
					
					
				}
				
				
				for (MismatchCounter c : useCounters)
					c.count(rd.getTotalCountsForConditions(null,mode));
				
				for (String object : this.objectCounter.objectsToCount.apply(readinfo))
					objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).count(rd.getTotalCountsForConditions(null,mode));
				
				if (this.intronCounter!=null && this.intronCounter.predicate.test(readinfo))
					for (String object : this.intronCounter.objectsToCount.apply(readinfo))
						objectCounter.computeIfAbsent(object, x->new MismatchCounter(object, cond, null, null)).count(rd.getTotalCountsForConditions(null,mode));
			}
			
		}
		
		synchronized (this.counter) {
			for (int i=0; i<counter.length; i++)
				this.counter.get(i).add(counter[i]);
		}
		
		if (numi!=null) 
			numiReport.accept(gene, numi.compute());
		
		return EI.wrap(objectCounter.values())
				.filter(g->g.reads.sum()>0)
				.map(r->r.toGeneData());
		
	}
	
	
	private boolean isConsistent(ImmutableReferenceGenomicRegion<Transcript> t,
			ImmutableReferenceGenomicRegion<AlignedReadsData> read) {
		
		// this is ugly, but gets the job done (leader sequence!):
//		if (t.getReference().getName().equals("NC_045512"))
//			return t.getRegion().intersects(read.getRegion());

		
		
		if (lenientOverlap) {
			int over = t.getRegion().intersect(read.getRegion()).getTotalLength();
			boolean leni = 2*over>=read.getRegion().getTotalLength() && t.getRegion().isIntronConsistent(read.getRegion());
			boolean strict = read.getData().isConsistentlyContained(read,t,0);
			if (strict && !leni)
				throw new RuntimeException(t+" "+read);
			return leni;
//			GenomicRegion tread = read.getRegion().extendBack(-3).extendFront(-3);
//			int exonlen = t.getRegion().intersect(tread).getTotalLength();
//			int intronlen = t.getRegion().invert().intersect(tread).getTotalLength();
//			return exonlen>intronlen*10;
		}
		
		if (t.getData().isCoding() && t.getData().get5Utr(t).getRegion().isEmpty() && t.getData().get3Utr(t).getRegion().isEmpty())
			return t.getRegion().intersects(read.getRegion()) && t.getRegion().isIntronConsistent(read.getRegion());
		return read.getData().isConsistentlyContained(read,t,0);
	}

	private boolean isBase(char c) {
		int ind = SequenceUtils.inv_nucleotides[c];
		return ind>=0 && ind<4;
	}


	private static class MismatchCounter {
		private String name;
		private HashMap<CounterKey,MutableDouble> counter = new HashMap<>();
		private AutoSparseDenseDoubleArrayCollector[] totalMM1;
		private AutoSparseDenseDoubleArrayCollector[] totalMM2;
		private AutoSparseDenseDoubleArrayCollector[] doublehits;
		private HashMap<IntPair, ReadData> binomData;
		private HashMap<IntPair, ReadData> binomOverlapData;
		private HashMap<IntPair, ReadData> binomBothData;
//		private HashMap<IntPair, ReadData> binomOutsideWhenDoubleData;

		private NumericArray reads;
		
		private Predicate<ReadInfo> predicate;
		private Function<ReadInfo,String[]> objectsToCount;
//		private boolean verbose = false;
		
		public MismatchCounter(MismatchCounter parent) {
			this.name = parent.name;
			this.predicate = parent.predicate;
			this.objectsToCount = parent.objectsToCount;
			this.reads = parent.reads.copy().clear();
			
			totalMM1 = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				totalMM1[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, parent.totalMM1[0].length()>>3),parent.totalMM1[0].length());
			totalMM2 = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				totalMM2[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, parent.totalMM1[0].length()>>3),parent.totalMM1[0].length());

			doublehits = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				doublehits[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, parent.totalMM1[0].length()>>3),parent.totalMM1[0].length());
			this.binomData = new HashMap<>();
			this.binomOverlapData = new HashMap<>();
			this.binomBothData = new HashMap<>();
//			this.binomOutsideWhenDoubleData = new HashMap<>();
		}
		

		public MismatchCounter(String name, int cond,
				Predicate<ReadInfo> predicate, Function<ReadInfo,String[]> objectsToCount) {
			this.name = name;
			reads = NumericArray.createMemory(cond, NumericArrayType.Double);
			
			totalMM1 = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				totalMM1[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond);
			totalMM2 = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				totalMM2[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond);
			
			doublehits = new AutoSparseDenseDoubleArrayCollector[16];
			for (int i=0; i<16; i++) 
				doublehits[i] = new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond);
			this.binomData = new HashMap<>();
			this.binomOverlapData = new HashMap<>();
			this.binomBothData = new HashMap<>();
//			this.binomOutsideWhenDoubleData = new HashMap<>();
			this.predicate = predicate;
			this.objectsToCount = objectsToCount;
		}
		
		public void count(NumericArray counts) {
			reads.add(counts);
		}
		
//		public MismatchCounter setVerbose(boolean verbose) {
//			this.verbose = verbose;
//			return this;
//		}

		public MismatchCounter createForThread() {
			return new MismatchCounter(this);
		}
		
		public void countBinomial(IntPair pp, AlignedReadsData rd, int d, ReadCountMode mode) {
			int ttotal = pp.a;
			int conv = pp.b;
			ReadData rrd = binomData.get(pp);
//			if (verbose){
//				System.out.println(name+" "+rrd);
//			}
			int cond = totalMM1[0].length();
			if (rrd==null) binomData.put(pp.cpy(), rrd=new ReadData(ttotal, conv, rd.addCountsForDistinct(d, new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond),mode)));
			else 
				rd.addCountsForDistinct(d, rrd.getCount(), mode);
//			if (verbose){
//				System.out.println(rrd);
//				System.out.println();
//			}
		}
		
		
		public void countBothBinomial(IntPair pp, AlignedReadsData rd, int d, ReadCountMode mode) {
			int ttotal = pp.a;
			int conv = pp.b;
			ReadData rrd = binomBothData.get(pp);
			int cond = totalMM1[0].length();
			if (rrd==null) binomBothData.put(pp.cpy(), new ReadData(ttotal, conv, rd.addCountsForDistinct(d, new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond),mode)));
			else 
				rd.addCountsForDistinct(d, rrd.getCount(), mode);
		}
		
		public void countDoubleHitBinomial(IntPair pp, AlignedReadsData rd, int d, ReadCountMode mode) {
			int ttotal = pp.a;
			int conv = pp.b;
			ReadData rrd = binomOverlapData.get(pp);
			int cond = totalMM1[0].length();
			if (rrd==null) binomOverlapData.put(pp.cpy(), new ReadData(ttotal, conv, rd.addCountsForDistinct(d, new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond),mode)));
			else 
				rd.addCountsForDistinct(d, rrd.getCount(), mode);
		}
		
//		public void countOutsideWhenDouble(IntPair pp, AlignedReadsData rd, int d, ReadCountMode mode) {
//			int ttotal = pp.a;
//			int conv = pp.b;
//			ReadData rrd = binomOutsideWhenDoubleData.get(pp);
//			int cond = totalMM1[0].length();
//			if (rrd==null) binomOutsideWhenDoubleData.put(pp.cpy(), new ReadData(ttotal, conv, rd.addCountsForDistinct(d, new AutoSparseDenseDoubleArrayCollector(Math.max(50, cond>>3),cond),mode)));
//			else 
//				rd.addCountsForDistinct(d, rrd.getCount(), mode);
//		}
		
		/**
		 * Opposite strand is for strand unspecific sequencing!
		 * @param indg
		 * @param indr
		 * @param rd
		 * @param d
		 * @param readpos
		 * @param overlap
		 * @param oppositeStrand
		 */
		public void countMismatch(int indg, int indr, AlignedReadsData rd, int d, int readpos, boolean firstread, boolean secondread, boolean oppositeStrand, ReadCountMode mode, boolean[] no4su) {
			if (!(firstread && secondread))
				rd.addCountsForDistinct(d, (secondread?totalMM2:totalMM1)[indg*4+indr], mode);
			CounterKey key = new CounterKey(indg, indr, readpos, firstread && secondread, oppositeStrand);
			MutableDouble arr = counter.get(key);
			if (arr==null) counter.put(key.cpy(), arr=new MutableDouble());
			arr.N+=rd.getTotalCountForDistinct(d, mode,(c,v)->no4su[c]?0:v);
		}
		
		/**
		 * Opposite strand is for strand unspecific sequencing!
		 * @param indg
		 * @param indr
		 * @param rd
		 * @param d
		 * @param readpos
		 * @param overlap
		 * @param oppositeStrand
		 */
		public void countTotal(int ind, AlignedReadsData rd, int d, int readpos, boolean firstread, boolean secondread, boolean oppositeStrand, ReadCountMode mode, boolean[] no4su) {
			if (!(firstread && secondread))
				rd.addCountsForDistinct(d, (secondread?totalMM2:totalMM1)[ind*4+ind], mode);
			
			CounterKey key = new CounterKey(ind, ind, readpos, firstread && secondread, oppositeStrand);
			MutableDouble arr = counter.get(key);
			if (arr==null) counter.put(key.cpy(), arr=new MutableDouble());
			arr.N+=rd.getTotalCountForDistinct(d, mode, (c,v)->no4su[c]?0:v);
		}
		
		public void countDoubleHit(int indg, int indr, AlignedReadsData rd, int d, ReadCountMode mode) {
			rd.addCountsForDistinct(d, doublehits[indg*4+indr], mode);
		}
		
		public void countDoubleTotal(int i, AlignedReadsData rd, int d, ReadCountMode mode) {
			rd.addCountsForDistinct(d, doublehits[i*4+i], mode);
		}
		
		public GeneData toGeneData() {
			return new GeneData(name,binomData.values().toArray(new ReadData[0]),binomOverlapData.values().toArray(new ReadData[0]),binomBothData.values().toArray(new ReadData[0]),reads);
		}
		
		/**
		 * Must be synchronized!
		 * @param thread
		 */
		public void add(MismatchCounter thread) {
			
			reads.add(thread.reads);
			
			for (int i=0; i<16; i++)
				doublehits[i].add(thread.doublehits[i]);
			for (int i=0; i<16; i++)
				totalMM1[i].add(thread.totalMM1[i]);
			for (int i=0; i<16; i++)
				totalMM2[i].add(thread.totalMM2[i]);
			
			for (CounterKey key : thread.counter.keySet()) {
				MutableDouble arr = counter.get(key);
				if (arr==null) counter.put(key.cpy(), arr=new MutableDouble());
				arr.N+=thread.counter.get(key).N;
			}
			
			for (IntPair p : thread.binomData.keySet()) {
				if (binomData.containsKey(p))
					binomData.get(p).add(thread.binomData.get(p));
				else
					binomData.put(p, thread.binomData.get(p));
			}
			
			for (IntPair p : thread.binomOverlapData.keySet()) {
				if (binomOverlapData.containsKey(p))
					binomOverlapData.get(p).add(thread.binomOverlapData.get(p));
				else
					binomOverlapData.put(p, thread.binomOverlapData.get(p));
			}
			for (IntPair p : thread.binomBothData.keySet()) {
				if (binomBothData.containsKey(p))
					binomBothData.get(p).add(thread.binomBothData.get(p));
				else
					binomBothData.put(p, thread.binomBothData.get(p));
			}
			
//			for (IntPair p : thread.binomOutsideWhenDoubleData.keySet()) {
//				if (binomOutsideWhenDoubleData.containsKey(p))
//					binomOutsideWhenDoubleData.get(p).add(thread.binomOutsideWhenDoubleData.get(p));
//				else
//					binomOutsideWhenDoubleData.put(p, thread.binomOutsideWhenDoubleData.get(p));
//			}
		}
	}
	
	public static class CounterKey {
		int genomic;
		int read;
		int readPos;
		boolean overlap;
		boolean readOppositeStrand;
		public CounterKey(int genomic, int read, int readPos, boolean overlap, boolean readOppositeStrand) {
			super();
			this.genomic = genomic;
			this.read = read;
			this.readPos = readPos;
			this.overlap = overlap;
			this.readOppositeStrand = readOppositeStrand;
		}
		public CounterKey cpy() {
			return new CounterKey(genomic, read, readPos, overlap, readOppositeStrand);
		}
		public char getGenomic() {
			return SequenceUtils.nucleotides[genomic];
		}
		public int getRead() {
			return SequenceUtils.nucleotides[read];
		}
		public int getReadPos() {
			return readPos;
		}
		public boolean isOverlap() {
			return overlap;
		}
		public boolean isReadOppositeStrand() {
			return readOppositeStrand;
		}
		@Override
		public String toString() {
			return "CounterKey [genomic=" + genomic + ", read=" + read + ", readPos=" + readPos + ", overlap=" + overlap
					+ ", readOppositeStrand=" + readOppositeStrand + "]";
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + genomic;
			result = prime * result + (overlap ? 1231 : 1237);
			result = prime * result + read;
			result = prime * result + (readOppositeStrand ? 1231 : 1237);
			result = prime * result + readPos;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			CounterKey other = (CounterKey) obj;
			if (genomic != other.genomic)
				return false;
			if (overlap != other.overlap)
				return false;
			if (read != other.read)
				return false;
			if (readOppositeStrand != other.readOppositeStrand)
				return false;
			if (readPos != other.readPos)
				return false;
			return true;
		}
		public boolean isTotal() {
			return genomic==read;
		}
		public CounterKey getTotalKey() {
			return new CounterKey(genomic, genomic, readPos, overlap, readOppositeStrand);
		}
		
		
	}
	
	private static class ReadInfo {
		private ImmutableReferenceGenomicRegion<String> gene;
		private ArrayList<ImmutableReferenceGenomicRegion<Transcript>> tr;
		
		private ImmutableReferenceGenomicRegion<AlignedReadsData> read;
		private String[] genes = nil;
		private String[] overlapgenes = nil;
		private String[] locations = nil;
		public ReadInfo(ImmutableReferenceGenomicRegion<String> gene,
				ArrayList<ImmutableReferenceGenomicRegion<Transcript>> tr) {
			this.gene = gene;
			this.tr = tr;
		}
		
		public boolean isOppositeStrand() {
			return !gene.getReference().getStrand().equals(read.getReference().getStrand());
		}
		
	}
	
	/** 
	 * the old rd is not usable anymore after this!
	 * @param rd
	 * @return
	 */
	public static ReadData[] collapse(ReadData[] rd) {
		HashMap<IntPair,AutoSparseDenseDoubleArrayCollector> c = new HashMap<>();
		IntPair key = new IntPair(0,0);
		for (ReadData r : rd) {
			AutoSparseDenseDoubleArrayCollector na = c.get(key.set(r.getConversions(), r.getTotal()));
			if (na==null) c.put(key.cpy(), na=r.getCount());
			else na.add(r.getCount());
		}
		return EI.wrap(c.keySet()).map(ip->new ReadData(ip.b, ip.a, c.get(ip))).toArray(ReadData.class);
	}
	
	private static class IntPair {
		int a;
		int b;
		public IntPair(int a, int b) {
			super();
			this.a = a;
			this.b = b;
		}
		public IntPair set(int a, int b) {
			this.a = a;
			this.b = b;
			return this;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + a;
			result = prime * result + b;
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			IntPair other = (IntPair) obj;
			if (a != other.a)
				return false;
			if (b != other.b)
				return false;
			return true;
		}
		public IntPair cpy(){
			return new IntPair(a, b);
		}
		
	}
	
	
	private MismatchCounter find(String name) {
		return EI.wrap(counter).filter(m->m.name.equals(name)).getUniqueResult(true, true);
	}
	
	public GeneData getExtData(String name) {
		return find(name).toGeneData();
	}
	
	public Collection<ReadData> getBinomData(String type) {
		if (type.length()>0)
			return EI.wrap(counter).filter(r->r.name.equals(type)).getUniqueResult(true, true).binomData.values();
		
		if (this.useAllReadsForModel)
			return EI.wrap(counter).filter(r->r.name.equals("All")).getUniqueResult(true, true).binomData.values();
		return objectCounter.binomData.values();
	}
	public Collection<ReadData> getBinomOverlapData(String type) {
		if (type.length()>0)
			return EI.wrap(counter).filter(r->r.name.equals(type)).getUniqueResult(true, true).binomOverlapData.values();
		
		if (this.useAllReadsForModel)
			return EI.wrap(counter).filter(r->r.name.equals("All")).getUniqueResult(true, true).binomOverlapData.values();
		return objectCounter.binomOverlapData.values();
	}
//	public Collection<ReadData> getOutsideOfOverlapWithDoubleData() {
//		return objectCounter.binomOutsideWhenDoubleData.values();
//	}

	public double getMismatches(int cond, char genomic, char read, boolean secondread) {
		if (locations==null) return getMismatches("Exonic",cond, genomic, read,secondread);
		return getMismatches("Locations",cond, genomic, read,secondread);
	}
	
	public List<String> getMismatchCategories() {
		return EI.wrap(counter).filter(mc->!mc.binomData.isEmpty() || !mc.binomOverlapData.isEmpty()).map(mc->mc.name).list();
	}
	
	public double getMismatches(String name, int cond, char genomic, char read, boolean secondread) {
		int indg = SequenceUtils.inv_nucleotides[genomic];
		int indr = SequenceUtils.inv_nucleotides[read];
		MismatchCounter mc = find(name);
		if (Math.min(indg, indr)>=0 && Math.max(indg, indr)<4)
			return (secondread?mc.totalMM2:mc.totalMM1)[indg*4+indr].get(cond);
		return Double.NaN;
	}
	
	public double getMismatches(String name, CounterKey k) {
		MismatchCounter mc = find(name);
		MutableDouble arr = mc.counter.get(k);
		if (arr==null) 
			return 0;
		return arr.N;
		
	}

	public double getCoverage(String name, CounterKey k) {
		MismatchCounter mc = find(name);
		MutableDouble total = mc.counter.get(k.getTotalKey());
		if (total==null) 
			return 0;
		return total.N;
		
	}
	
	public double getDoubleHits(String name, int cond, char genomic, char read) {
		int indg = SequenceUtils.inv_nucleotides[genomic];
		int indr = SequenceUtils.inv_nucleotides[read];
		MismatchCounter mc = find(name);
		if (Math.min(indg, indr)>=0 && Math.max(indg, indr)<4)
			return mc.doublehits[indg*4+indr].get(cond);
		return Double.NaN;
	}
	public double getDoubleHitCoverage(String name, int cond, char genomic) {
		int indg = SequenceUtils.inv_nucleotides[genomic];
		MismatchCounter mc = find(name);
		if (indg>=0 && indg<4)
			return mc.doublehits[indg*4+indg].get(cond);
		return Double.NaN;
	}
	
	
	public double getCoverage(String name, int cond, char genomic, boolean secondread) {
		int indg = SequenceUtils.inv_nucleotides[genomic];
		MismatchCounter mc = find(name);
		if (indg>=0 && indg<4)
			return (secondread?mc.totalMM2:mc.totalMM1)[indg*4+indg].get(cond);
		return Double.NaN;
	}
	
	public ExtendedIterator<CounterKey> getMismatchPositionKeys(String name) {
		MismatchCounter mc = find(name);
		return EI.wrap(mc.counter.keySet()).filter(c->!c.isTotal());
	}
	
	
//	public Collection<ReadData> getIntronData() {
//		return cintronData.values();
//	}
//
//	public double getIntronicConversionRate(int cond) {
//		double conv = 0;
//		double total = 0;
//		for (ReadData rd : getIntronData()) {
//			conv+=rd.getConversions()*rd.getCount().getDouble(cond);
//			total+=rd.getTotal()*rd.getCount().getDouble(cond);
//		}
//		return conv/total;
//	}
	
	private static class MutableSingletonGenomicRegion implements GenomicRegion {
		
		int p;


		@Override
		public int getNumParts() {
			return 1;
		}

		@Override
		public int getStart(int part) {
			return p;
		}

		@Override
		public int getEnd(int part) {
			return p+1;
		}
		
		@Override
		public String toString() {
			return toRegionString();
		}
		
		
	}

	

	
}
