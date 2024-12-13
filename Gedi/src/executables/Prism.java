package executables;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLongArray;
import java.util.function.BiFunction;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToIntFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import cern.colt.bitvector.BitVector;
import gedi.app.Gedi;
import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.annotation.NameAnnotation;
import gedi.core.data.annotation.ScoreNameAnnotation;
import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.reference.Strandness;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.intervalTree.MemoryIntervalTreeStorage;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.FunctorUtils.FilteredIterator;
import gedi.util.LogUtils.LogMode;
import gedi.util.ReflectionUtils;
import gedi.util.RunUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.charsequence.CharDag;
import gedi.util.datastructure.charsequence.CharIterator;
import gedi.util.datastructure.charsequence.TranslateCharDagVisitor;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.datastructure.tree.Trie.AhoCorasickResult;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.functions.StringIterator;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.fasta.DefaultFastaHeaderParser;
import gedi.util.io.text.fasta.FastaEntry;
import gedi.util.io.text.fasta.FastaFile;
import gedi.util.io.text.fasta.FastaHeaderParser;
import gedi.util.io.text.tsv.formats.Bed;
import gedi.util.io.text.tsv.formats.BedEntry;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.math.stat.counting.Counter;
import gedi.util.mutable.MutableMonad;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.GenomicParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;
import gedi.util.r.RRunner;
import gedi.util.sequence.Alphabet;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.Progress;


@SuppressWarnings("unused")
public class Prism {
	
	

	private static String getChangelog() {
		return "1.1.1:\n"
				+ "introduced the -dmod parameter\n\n"
				+ "1.1.1a:\n"
				+ "fixed a minor bug with peptides matching into a longer variant (but not completely covering it)\n\n"
				+ "1.1.2:\n"
				+ "Prism now recognizes frameshifted peptdes downstream of an indel (if consistent) as CDS; only a single indel is considered!\n\n"
				+ "1.1.2a:\n"
				+ "Fixed bugs for the frameshifted peptides (when Stop codon is within insertion or downstream of annotated transcript)!\n\n"				
				+ "1.1.2b:\n"
				+ "Fixed another bug introduced with the frameshifted peptides!\n\n"
				+ "1.1.3:\n"
				+ "Added support for RNA-seq reads; Delta next output & fixed a minor bug (incorrect output if a peptide ends exactly with the end of a variant)!\n\n"
				+ "1.1.3a:\n"
				+ "Fixed exception with two fastq files\n\n"
				+ "1.1.3b:\n"
				+ "Fixed bug with annotating RNA-seq reads\n\n"
				+ "1.1.3c:\n"
				+ "Fixed another bug for RNA-seq reads\n\n"
				+ "1.1.3d:\n"
				+ "Fixed yet another bug for RNA-seq reads\nAdded 'Top location count (no decoy)' (which is the number of equally good choices for the location before prioritizing by reads and/or ORF length)\n\n"
				+ "1.1.4:\n"
				+ "Output of DNA sequences and fastq files for category READS; Reads can now be priority 1\n\n"
				+ "1.1.4a:\n"
				+ "Genome may now contain lowercase letters (e.g. T2T genome). Fixed issues with determining start positions for deletions\n\n"
				+ "1.1.4b:\n"
				+ "Fixed issues with lower case letters within variants\n\n"
				+ "1.1.4c:\n"
				+ "Fixed issues with data frame names in specific R versions\n\n"
				+ "1.1.4d:\n"
				+ "Fixed issue with multiple alleles in vcf; non-matching sequences are now only reported as warnings; : replaced by - in fastq output!\n\n"
				+ "1.1.4g:\n"
				+ "Exclude incomplete transcripts from annotating peptides.\n\n"
				+ "1.1.4h:\n"
				+ "Fixed bed output bug (missing netMHC score).\n\n"
				+ "1.1.4i:\n"
				+ "Fixed underlying bug for bed output (comma in frameshift ORF shifting table row).\n\n"
				+ "1.1.5:\n"
				+ "Fixed read 2 bug; added into intron categories.\n\n"
				+ "1.1.5a:\n"
				+ "Set lower default priority to OtherIntoIntron.\n\n"
				+ "1.1.6:\n"
				+ "Fixed prioritization of genomic categories (prior to that if the same location was e.g. CDS for one and 5UTR for another transcript, CDS always won (no matter what priorities were defined).\n\n"
				+ "1.1.7:\n"
				+ "If there are non-decoys in a lower priority category, they are now preferred over decoys in higher categories.\n\n";
	}


	public static void main(String[] args) throws IOException {
		
		Gedi.startup(true,LogMode.Normal,"Peptide-PRISM");
		
		FindGenomicPeptidesParameterSet params = new FindGenomicPeptidesParameterSet();
		GediProgram pipeline = EI.wrap(args).filter(p->p.equals("-anchor")).count()>0?
				GediProgram.create("Prism",
						new PeptidesCountProgram(params)
					)
				: GediProgram.create("Prism",
					new FindGenomicPeptidesProgram(params),
					new ToPeptideListProgram(params),
					new ComputeFdrProgram(params),
					new NetMHCProgram(params),
					new SrrCalcProgram(params),
					new AnnotateProgram(params),
					new BedOutProgram(params),
					new SearchReadsProgram(params)
				//	new PlotProgram(params)
	//				new ProcessUnidentifiedProgram(params)
				);
		pipeline.setChangelog(getChangelog());
		
		GediProgram.run(pipeline, params.paramFile, params.runtimeFile, new CommandLineHandler("Prism","Prism identifies a list of peptides in the genomic 6-frame translation.",args));

	}
	
	public static class PeptideOrf implements Comparable<PeptideOrf>{
		String dna;
		int[] startPrioPos;
		int[] startNextPos;
		GenomicRegion pepLocInDna;
		public PeptideOrf(String dna, GenomicRegion pepLocInDna, int[][] startPos) {
			this.dna = dna;
			this.startPrioPos = startPos==null?null:startPos[0];
			this.startNextPos = startPos==null?null:startPos[1];
			this.pepLocInDna = pepLocInDna;
		}
		
		@Override
		public int compareTo(PeptideOrf data) {
			if (startPrioPos==null && data.startPrioPos==null) return 0;
			if (startPrioPos==null) return 1;
			if (data.startPrioPos==null) return -1;
			int re = Integer.compare(startPrioPos[0], data.startPrioPos[0]);
			if (re==0)
				re = Integer.compare(data.startPrioPos[1], startPrioPos[1]);
			return re;
		}

		public String getTisPrio() {
			return startPrioPos==null?"-":dna.substring(startPrioPos[1], startPrioPos[1]+3);
		}
		
		public String getTis() {
			return startNextPos==null?"-":dna.substring(startNextPos[1], startNextPos[1]+3);
		}
		
		public String getSeqPrio() {
			return startPrioPos==null?null:dna.substring(startPrioPos[1]);
		}
		
		public String getSeq() {
			return startNextPos==null?null:dna.substring(startNextPos[1]);
		}

		public int getLength() {
			return startPrioPos==null?-1:dna.length()-startPrioPos[1];
		}

		public int getNterm() {
			return startPrioPos==null?-1:(pepLocInDna.getStart()-startPrioPos[1])/3;
		}

		public int getCterm() {
			return startPrioPos==null?-1:(dna.length()-pepLocInDna.getEnd())/3-1; // stop codon
		}
	}
	
	public static class PeptideToRegion implements Comparable<PeptideToRegion>{
		Genomic g;
		ReferenceGenomicRegion<Transcript> transcript;
		ImmutableReferenceGenomicRegion<Transcript> region;
		ImmutableReferenceGenomicRegion<?> pep;
		
		public PeptideToRegion(Genomic g,ImmutableReferenceGenomicRegion<Transcript> region,
					ImmutableReferenceGenomicRegion<?> pep) {
			this.g = g;

			this.transcript = region.getData()==null?null:g.getTranscriptMapping().apply(region.getData().getTranscriptId());
			this.region = region;
			this.pep = pep;
		}
		
		@Override
		public int compareTo(PeptideToRegion data) {
			boolean tin = region.contains(pep);
			boolean din = data.region.contains(data.pep);
			
			if (tin!=din)
				return tin?-1:1;
			
			boolean tccds = region.getData()!=null && g.getTranscriptTable("source") !=null && "CCDS".equals(g.getTranscriptTable("source").apply(region.getData().getTranscriptId()));
			boolean dccds = data.region.getData()!=null && g.getTranscriptTable("source") !=null && "CCDS".equals(g.getTranscriptTable("source").apply(data.region.getData().getTranscriptId()));
			
			if (tccds!=dccds)
				return tccds?-1:1;
			
			return Integer.compare(data.region.getRegion().getTotalLength(), region.getRegion().getTotalLength());
		}
		
		public ImmutableReferenceGenomicRegion<Transcript> getRegion() {
			return region;
		}

		public int getNterm() {
			int s = pep.map(0);
			if (transcript!=null) {
				int intrans = transcript.induce(s);
				ArrayGenomicRegion regInTrans = transcript.induce(region.getRegion());
				return intrans-regInTrans.getStart();
			}
			return region.induceMaybeOutside(s);
		}

		public int getCterm() {
			int s = pep.map(pep.getRegion().getTotalLength()-1);
			if (transcript!=null) {
				int intrans = transcript.induce(s);
				ArrayGenomicRegion regInTrans = transcript.induce(region.getRegion());
				return regInTrans.getStop()-intrans;
			}
			return region.getRegion().getTotalLength()-1-region.induceMaybeOutside(s);
		}
		
		public int getIntronIndex() {
			for (ImmutableReferenceGenomicRegion<Transcript> t : g.getGenes().ei(region)
													.unfold(gene->g.getTranscripts().ei(gene)).loop()) {
				
				for (int i=0; i<t.getRegion().getNumParts()-1; i++)
					if (region.getRegion().equals(t.getRegion().getIntron(i)))
						return t.getReference().isPlus()?i:(t.getRegion().getNumParts()-2)-i;
			}
			return -1;
		}
		
		public String getUpstream(int l) {
			return SequenceUtils.translate(g.getSequenceSave(pep.getUpstream(l)));
		}
		
		public String getDownstream(int l) {
			return SequenceUtils.translate(g.getSequenceSave(pep.getDownstream(l)));
		}

	}
	
	
	public static MemoryIntervalTreeStorage<Transcript> intoIntron = null;
	public static int[] discoverIntoIntron(Genomic genomic, Progress progress) {
		int[] re = {0,0};
		intoIntron = new MemoryIntervalTreeStorage<Transcript>(Transcript.class);
		
		Trie<Void> stopCodons = new Trie<Void>();
		stopCodons.put("TAA", null);
		stopCodons.put("TGA", null);
		stopCodons.put("TAG", null);
		BitVector frameFound = new BitVector(3);
		
		for (ImmutableReferenceGenomicRegion<Transcript> tr : genomic.getTranscripts().ei().progress(progress, (int)genomic.getTranscripts().size(), tr->tr.toLocationString()).loop()) {
			GenomicRegion cds = tr.getData().isCoding()?tr.getData().getCds(tr.getReference(), tr.getRegion()):new ArrayGenomicRegion();
			for (int ex=0; ex<tr.getRegion().getNumParts()-1; ex++) {
				GenomicRegion upstream = tr.getReference().isMinus()?
						tr.getRegion().intersect(tr.getRegion().getStart(ex+1), tr.getRegion().getEnd()):
							tr.getRegion().intersect(tr.getRegion().getStart(), tr.getRegion().getEnd(ex));
				GenomicRegion upstreamCds = upstream.intersect(cds);
				
				GenomicRegion intron = tr.getRegion().getIntron(ex);
				int codingFrame = upstreamCds.isEmpty()?-1:(3-upstreamCds.getTotalLength()%3)%3;
				
				CharSequence intronSeq = genomic.getSequence(tr.getReference(), intron);
				for (AhoCorasickResult<Void> res : stopCodons.iterateAhoCorasick(intronSeq).loop()) {
					int frame = res.getStart()%3;
					if (!frameFound.getQuick(frame)) {
						// this is the first stop codon in this frame!
						frameFound.putQuick(frame, true);
						
						ArrayGenomicRegion intronret = upstream.union(tr.getReference().isMinus()?
								new ArrayGenomicRegion(upstream.getStart()-res.getStart(),upstream.getStart()):
									new ArrayGenomicRegion(upstream.getEnd(),upstream.getEnd()+res.getStart()));
						Transcript trans;
						if (codingFrame==frame) {
							re[0]+=res.getStart();
							if (tr.getReference().isMinus())
								trans = new Transcript(tr.getData().getGeneId(), tr.getData().getTranscriptId(), intronret.getStart(), tr.getData().getCodingEnd());
							else
								trans = new Transcript(tr.getData().getGeneId(), tr.getData().getTranscriptId(), tr.getData().getCodingStart(), intronret.getEnd());
						} else {
							re[1]+=res.getStart();
							trans = new Transcript(tr.getData().getGeneId(), tr.getData().getTranscriptId(), -1, -1);
						}
						intoIntron.add(tr.getReference(), intronret, trans);
						if (frameFound.cardinality()==3) break;
					}
						
				}
				frameFound.clear();
			}
		}
		
//		try {
//			new CenteredDiskIntervalTreeStorage<>("test.cit", Transcript.class).fill(intoIntron);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
		
		return re;
	}
	
	
	
	public static MemoryIntervalTreeStorage<Transcript> buildFrameshiftCds(Genomic genomic, GenomicRegionStorage<NameAnnotation> vars) {
		MemoryIntervalTreeStorage<String[]> indels = new MemoryIntervalTreeStorage<>(String[].class);
		for (ReferenceGenomicRegion<NameAnnotation> vv : vars.ei().loop()) {
			String[] fromto = StringUtils.split(vv.getData().getName(), '>');
			if (fromto[1].contains(",")) fromto[1] = fromto[1].substring(0,fromto[1].indexOf(',')); // if multi-var, only use the first one!
			if ((fromto[1].length()-fromto[0].length())%3!=0) {
				indels.add(vv.getReference().toPlusStrand(), vv.getRegion(), fromto);
				indels.add(vv.getReference().toMinusStrand(), vv.getRegion(), fromto);
			}
		}
		
		MemoryIntervalTreeStorage<Transcript> re = new MemoryIntervalTreeStorage<>(Transcript.class);

		for (ReferenceGenomicRegion<Transcript> t : genomic.getTranscripts().ei().
				filter(e->e.getData().isCoding()).loop()) {
//			if (t.getData().getGeneId().equals("ENSG00000147889"))
//				indels.ei(t).print();
			
			for (ReferenceGenomicRegion<String[]> vv : indels.ei(t).filter(vv->t.getData().getCds(t).getRegion().containsUnspliced(vv.getRegion())).loop()) {
				
				GenomicRegion cdsutr = t.map(new ArrayGenomicRegion(t.getData().get5UtrLength(t.getReference(), t.getRegion()),t.getRegion().getTotalLength()));
				
				int chrlen = genomic.getLength(t.getReference().getName());
				if (t.getReference().isPlus()) 
					cdsutr=cdsutr.extendBack(Math.min(chrlen-cdsutr.getEnd(), 3000));
				else 
					cdsutr=cdsutr.extendFront(Math.min(cdsutr.getStart(), 3000));
				// in case there is no stop codon add 3000 to the back!
				
				String seq = genomic.getSequence(vv.getReference().toPlusStrand(), cdsutr).toString().toUpperCase();
				int varStart = cdsutr.induce(vv.getRegion().getStart());
				int varEnd = cdsutr.induce(vv.getRegion().getStop())+1;
				// this is in genome position space, not 5'->3'!
				
				if (!seq.substring(varStart, varEnd).equals(vv.getData()[0]))
					throw new RuntimeException("Variant sequence does not match genome: "+vv.getReference()+":"+vv.getRegion()+" "+Arrays.toString(vv.getData())+" Transcript: "+t);
				
				String replaced = seq.substring(0,varStart)+vv.getData()[1]+seq.substring(varEnd);
				if (vv.getReference().isMinus()) replaced=SequenceUtils.getDnaReverseComplement(replaced);
				int upstream = vv.getReference().isPlus()?varStart:seq.length()-varEnd;
				int stop = SequenceUtils.translate(replaced).indexOf('*')*3;
				if (stop==-1) stop=cdsutr.getTotalLength()-3;
				else if (stop<=upstream+vv.getData()[1].length()) // STOP inside of insertion
					stop=upstream;
				else stop=stop-vv.getData()[1].length()+vv.getData()[0].length();
				
				//obtain only the part after the frameshift!
				GenomicRegion newCds;
				if (stop==upstream)
					newCds = new ArrayGenomicRegion();
				else if (t.getReference().isPlus()) {
					int downstreamOfIndel = Math.max(0,stop-varEnd); // if stop within insertion!
					downstreamOfIndel = downstreamOfIndel/3*3;
					int startOfNewCds = stop-downstreamOfIndel;
					
					newCds = cdsutr.map(startOfNewCds, stop+3);
				}
				else {
					int tstop = cdsutr.getTotalLength()-(stop+3);
					int downstreamOfIndel = Math.max(0,varStart-tstop);
					downstreamOfIndel = downstreamOfIndel/3*3;
					
					newCds = cdsutr.map(tstop,tstop+downstreamOfIndel);
				}
				
				if (newCds.getTotalLength()==0) {
					// either: var goes over indel (then this is handeled as direct frameshift peptide)
					// or: the ORF is invalid
					if (!new ImmutableReferenceGenomicRegion<>(t.getReference(), cdsutr).map(new ArrayGenomicRegion(stop,stop+3)).intersects(vv.getRegion()) && SequenceUtils.checkCompleteCodingTranscript(genomic, t, 0, 0, true)) 
						throw new RuntimeException("Cannot be! "+vv.toString());
					continue;
				}
				
				
//				if (!t.getRegion().containsUnspliced(newCds)) 
//					throw new RuntimeException("Cannot be!"+vv.toString());
				// can be, if the next stop codon is downstream of the gene
				
				
				String suff = "@"+vv.toLocationString()+"/"+vv.getData()[0]+">"+vv.getData()[1];
//				ImmutableReferenceGenomicRegion<Transcript> ref = new ImmutableReferenceGenomicRegion<>(t.getReference(),t.getRegion(),new Transcript(t.getData().getGeneId()+suff, t.getData().getTranscriptId()+suff, newCds.getStart(), newCds.getEnd()));
				ImmutableReferenceGenomicRegion<Transcript> ref = new ImmutableReferenceGenomicRegion<>(t.getReference(),t.getReference().isPlus()?t.getRegion().extendBack(3000):t.getRegion().extendFront(3000),new Transcript(t.getData().getGeneId()+suff, t.getData().getTranscriptId()+suff, newCds.getStart(), newCds.getEnd()));
				if (!ref.getData().getCds(ref).getRegion().equals(newCds))
					throw new RuntimeException("Cannot be!");
				re.add(ref);
//				if (t.getData().getGeneId().equals("ENSG00000147889"))
//					System.out.println(ref.getData().getCds(ref).toLocationString()+" "+ref.getData().getCds(ref).induce(ImmutableReferenceGenomicRegion.parse("9-:21971000-21971027").getRegion()));
			}
		}
		return re;
	}
	
	
	public static List<ImmutableReferenceGenomicRegion<PeptideOrf>> getOrfs(ImmutableReferenceGenomicRegion<?> pep, Genomic g, Trie<Integer> starts, boolean shortest) {
		
		if (pep.getReference().getName().startsWith("REV_"))
			pep = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(pep.getReference().getName().substring(4),pep.getReference().getStrand()), pep.getRegion());
		
		ImmutableReferenceGenomicRegion<?> upep = pep;
		
		ArrayList<ImmutableReferenceGenomicRegion<PeptideOrf>> re = new ArrayList<>();
		if (pep.getRegion()==null || pep.getRegion().getTotalLength()%3!=0) return re;
		int ex = 99;
		int lu = -1;
		int ld = -1;
		
		do  {
			int uex = ex;
			if (pep.getRegion().extendFront(uex).extendBack(uex).subtract(new ArrayGenomicRegion(0,g.getLength(upep.getReference().getName()))).getTotalLength()>0)
				return re;
			
			String seq = SequenceUtils.translate(g.getSequence(pep.toMutable().transformRegion(r->r.extendFront(uex).extendBack(uex))));
			
			int up = seq.substring(0, ex/3).lastIndexOf('*');
			ld = seq.indexOf('*', ex/3)*3+3-ex-pep.getRegion().getTotalLength();
			lu=up<0?-1:(ex/3-up)*3-3;
			ex=ex*3;
		} while ((lu<0 || ld<0) && pep.getRegion().getStart()-ex/3>=0 && pep.getRegion().getEnd()+ex/3<g.getLength(pep.getReference().getName()));
		
		if (lu>=0 && ld>=0) {
			MutableReferenceGenomicRegion<?> orf = pep.toMutable().extendRegion(lu,ld);
			String seq = g.getSequence(orf).toString();
			
			re.add(new ImmutableReferenceGenomicRegion<>(orf.getReference(), orf.getRegion(), new PeptideOrf(seq,new ArrayGenomicRegion(lu,orf.getRegion().getTotalLength()-ld),findStart(seq,lu,starts, shortest))));
		}
		
		for (ImmutableReferenceGenomicRegion<Transcript> t : g.getTranscripts().ei(pep).filter(t->t.getRegion().containsUnspliced(upep.getRegion())).loop()) {
			int pepPosAA = t.induce(pep.getRegion()).getStart()/3;
			int offset = t.induce(pep.getRegion()).getStart()%3;
			
			String dna = g.getSequence(t).toString().substring(offset);
			String seq = SequenceUtils.translate(dna);
			
			int up = seq.substring(0, pepPosAA).lastIndexOf('*');
			if (up==-1) {
				if (t.getData().isCoding() && t.getData().getCds(t).getRegion().containsUnspliced(pep.getRegion()) && t.getData().getCds(t).induce(pep.getRegion()).getStart()%3==0)
					up = -1;
				else
					up = 4;
			}
			lu=(pepPosAA-up)*3-3;
			ld = seq.indexOf('*', pepPosAA)*3+3-pepPosAA*3-pep.getRegion().getTotalLength();
			
			if (lu>=0 && ld>=0) {
				ArrayGenomicRegion xorf = new ArrayGenomicRegion(up*3+3,pepPosAA*3+pep.getRegion().getTotalLength()+ld);
				seq = xorf.extractSequence(dna);
				
				re.add(new ImmutableReferenceGenomicRegion<>(t.getReference(), t.toMutable().extendRegion(-offset, 0).map(xorf), new PeptideOrf(seq,new ArrayGenomicRegion(lu,xorf.getTotalLength()-ld),findStart(seq,lu,starts, shortest))));
			}
		}
		
		return re;
	}
	
	/**
	 * [start codon prio,start pos] [start codon next,start pos]
	 * @param seq
	 * @param pep
	 * @param starts
	 * @return
	 */
	private static int[][] findStart(String seq, int pep, Trie<Integer> starts, boolean shortest) {
		ArrayList<int[]> list = starts.iterateAhoCorasick(seq)
				.filter(res->res.getStart()%3==0 & res.getStart()<=pep)
				.map(res->new int[]{res.getValue(),res.getStart()})
				.list();
		if (list.size()==0) return null;
		list.sort((a,b)->{
			int re = Integer.compare(a[0], b[0]); // first by start codon priority
			if (re==0)
				re = shortest?Integer.compare(b[1], a[1]):Integer.compare(a[1], b[1]); // then by shorter ORF
			return re;
		});
		int[] re1 = list.get(0);
		if (shortest)
			list.sort((a,b)->Integer.compare(b[1], a[1]));
		else
			list.sort((a,b)->Integer.compare(a[1], b[1]));
		return new int[][] {re1,list.get(0)};
	}
	
	public static MutablePair<Category,List<ImmutableReferenceGenomicRegion<Transcript>>> annotateLocation(Genomic g, AnnotationEnum aenum, ImmutableReferenceGenomicRegion<String> ahit, MemoryIntervalTreeStorage<Transcript> extraCds) {
		if (ahit.getReference().equals(readDummyRegion.getReference()) || ahit.getReference().equals(readDecoyDummyRegion.getReference())) return new MutablePair<>(Category.READS,null);
		if (g==null) return new MutablePair<>(Category.CDS,null);
		
		if (ahit.getReference().getName().startsWith("REV_"))
			ahit = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(StringUtils.removeHeader(ahit.getReference().getName(), "REV_"),ahit.getReference().getStrand()), ahit.getRegion(),ahit.getData());
		
		ImmutableReferenceGenomicRegion<String> hit = ahit;
		ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs = g.getTranscripts().ei(hit).list();
		MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf = new MutableMonad<>(new ArrayList<>());
		
		
		if (ahit.getData().endsWith("RNA-seq"))
			return new MutablePair<>(Category.RNASEQ,null);
	
		
		if (ahit.getData().endsWith("Extra"))
			return new MutablePair<>(Category.Extra,null);
	
		int sind = ahit.getData().indexOf('>');
		int cind = ahit.getData().indexOf(':');
		if (cind>0 && cind+2==sind)
			return new MutablePair<>(Category.Substitution,null);
		
		if (EI.wrap(trs)
				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
				.map(tr->tr.getData().getCds(tr))
				.sideEffect(r->buf.Item.add(r.toImmutable()))
				.count()>0
				&& hit.getRegion().getTotalLength()%3!=0)
				return new MutablePair<>(Category.Frameshift,buf.Item);
		
		if (StringUtils.countChar(ahit.getData(),':')>1)
			return new MutablePair<>(Category.PeptideSpliced,null);
	
	
		for (Annotation a : aenum.values) {
			for (Category c : a.categories) {
				
				if (c.matches(g, hit, extraCds, trs, buf)) {
					buf.Item.clear();
					return new MutablePair<Category, List<ImmutableReferenceGenomicRegion<Transcript>>>(c, buf.Item);
				}
				
			}
		}
		
//		buf.Item.clear();
//		if (EI.wrap(trs)
//			.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
//			.map(tr->tr.getData().getCds(tr))
//			.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
//			.sideEffect(r->buf.Item.add(r.toImmutable()))
//			.count()>0)
//			return new MutablePair<>(Category.CDS,buf.Item);
//
//		buf.Item.clear();
//		if (extraCds!=null && extraCds.ei(hit)
//			.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
//			.map(tr->tr.getData().getCds(tr))
//			.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
//			.sideEffect(r->buf.Item.add(r.toImmutable()))
//			.count()>0)
//			return new MutablePair<>(Category.CDS,buf.Item);

//		buf.Item.clear();
//		if (EI.wrap(trs)
//				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
//				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
//				.map(tr->tr.getData().get5Utr(tr))
//				.filter(tr->tr.getRegion().intersects(hit.getRegion()) && tr.getRegion().isIntronConsistent(hit.getRegion()))
//				.sideEffect(r->buf.Item.add(r.toImmutable()))
//				.count()>0)
//				return new MutablePair<>(Category.UTR5,buf.Item);
		
//		buf.Item.clear();
//		if (EI.wrap(trs)
//				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
//				.map(tr->tr.getData().getCds(tr))
//				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
//				.sideEffect(r->buf.Item.add(r.toImmutable()))
//				.count()>0)
//				return new MutablePair<>(Category.OffFrame,buf.Item);

//		buf.Item.clear();
//		if (EI.wrap(trs)
//				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
//				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
//				.map(tr->tr.getData().get3Utr(tr))
//				.filter(tr->tr.getRegion().intersects(hit.getRegion()))
//				.sideEffect(r->buf.Item.add(r.toImmutable()))
//				.count()>0)
//				return new MutablePair<>(Category.UTR3,buf.Item);
		
		
		
//		buf.Item.clear();
//		if (EI.wrap(trs)
//				.filter(tr->!tr.getData().isCoding() || !SequenceUtils.checkCompleteCodingTranscript(g, tr))
//				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
//				.sideEffect(r->buf.Item.add(r.toImmutable()))
//				.count()>0)
//				return new MutablePair<>(Category.ncRNA,buf.Item);
		
		if (EI.wrap(trs)
				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
				.map(tr->tr.getData().getCds(tr))
				.filter(tr->tr.getRegion().contains(hit.getRegion()))
				.count()>0)
			throw new RuntimeException(hit.toLocationString());
		
		
		if (EI.wrap(trs)
				.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
				.count()>0)
				throw new RuntimeException(hit.toLocationString());
		
//		buf.Item.clear();
//		if (intoIntron!=null && EI.wrap(intoIntron.ei(hit))
//			.filter(tr->tr.getData().isCoding())
//			.map(tr->tr.getData().getCds(tr))
//			.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
//			.sideEffect(r->buf.Item.add(new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), r.getRegion())))
//			.count()>0)
//			return new MutablePair<>(Category.CDSintoIntron,buf.Item);

//		buf.Item.clear();
//		if (intoIntron!=null && EI.wrap(intoIntron.ei(hit))
//				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
//				.sideEffect(r->buf.Item.add(new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), r.getRegion())))
//				.count()>0)
//				return new MutablePair<>(Category.OtherIntoIntron,buf.Item);

		
		
//		buf.Item.clear();
//		if (g.getGenes().ei(hit)
//				.filter(gene->gene.contains(hit))
//				.count()>0) {
//			
//			g.getGenes().ei(hit)
//				.unfold(gene->g.getTranscripts().ei(gene))
//				.unfold(t->EI.wrap(t.getRegion().invert().getParts()))
//				.map(p->p.asRegion())
//				.filter(p->p.intersects(hit.getRegion()))
//				.map(p->new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), p))
//				.toCollection(buf.Item);
//			
//				return new MutablePair<>(Category.Intronic,buf.Item);
//		}
		
//		buf.Item.clear();
//		ArrayList<GenomicRegion> left = g.getTranscripts().getTree(hit.getReference()).getIntervalsLeftNeighbor(hit.getRegion().getStart(), hit.getRegion().getStart(), new ArrayList<>());
//		ArrayList<GenomicRegion> right = g.getTranscripts().getTree(hit.getReference()).getIntervalsRightNeighbor(hit.getRegion().getStop(), hit.getRegion().getStop(), new ArrayList<>());
//		if (!left.isEmpty() && !right.isEmpty()) {
//			buf.Item.add(new ImmutableReferenceGenomicRegion<>(hit.getReference(), new ArrayGenomicRegion(left.get(0).getEnd(),right.get(0).getStart())));
//		}
//		
//		return new MutablePair<>(Category.Intergenic,buf.Item);
	
		
		buf.Item.clear();
		return new MutablePair<>(Category.Unknown,buf.Item);
		
	}
	
	
	
	public static String getMutationSequence(Genomic g, String chr, int s, int e, String ref,
			String alt, int flank) {
		ImmutableReferenceGenomicRegion<Void> p = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(chr, true),new ArrayGenomicRegion(s-1,e));
		ImmutableReferenceGenomicRegion<Void> m = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(chr, false),new ArrayGenomicRegion(s-1,e));
		
		return g.getTranscripts().ei(p).chain(g.getTranscripts().ei(m))
			.filter(tr->tr.getData().isCoding())
			.map(tr->tr.getData().getCds(tr))
			.filter(cds->cds.getRegion().contains(p.getRegion()))
			.map(cds->getFlanking(g,cds, p, ref, alt, cds.getReference().isPlus(),flank)).first();
	}
	
	private static final int getScore(HeaderLine h, String[] f) {
		int alc = h.get("ALC (%)",-1);
		if (alc!=-1) return Integer.parseInt(f[alc]);
		
		int lgp = h.get("-10lgP",-1);
		if (lgp!=-1) return (int)(100*Double.parseDouble(f[lgp]));
		
		return 0;
	}
	
	
	private static String getFlanking(Genomic g, MutableReferenceGenomicRegion<?> ref,
			ImmutableReferenceGenomicRegion<?> pos, String from, String to, boolean plus, int flank) {
		
		int p = ref.induce(pos.getRegion()).getStart();
		if (!plus)
			to = SequenceUtils.getDnaComplement(to);
		
		String dna = g.getSequence(ref).toString();
		String mut = dna.substring(0,p)+to+dna.substring(p+1);
		String aa = StringUtils.repeat(' ', -Math.min(0, p/3-flank))
				+SequenceUtils.translate(mut).substring(Math.max(0, p/3-flank),Math.min(mut.length()/3, p/3+flank))
				+StringUtils.repeat(' ', Math.max(mut.length()/3, p/3+flank)-mut.length()/3);
		return aa;
	}

	public static class Affinity implements Comparable<Affinity> {
		public String allele;
		public String sequence;
		public double perrank;
		
		public Affinity(String allele, String sequence, double perrank) {
			this.allele = allele;
			this.sequence = sequence;
			this.perrank = perrank;
		}

		@Override
		public int compareTo(Affinity o) {
			return Double.compare(perrank, o.perrank);
		}
		
		@Override
		public String toString() {
			return sequence+" "+allele+" "+perrank;
		}
	}

	public static final Category[] defaultCategories = {
			Category.CDS,
			Category.RNASEQ,
			Category.UTR5,
			Category.OffFrame,
			Category.CDSintoIntron,
			Category.UTR3,
			Category.ncRNA,
			Category.OtherIntoIntron,
			Category.Extra,
			Category.READS,
			Category.Intronic,
			Category.Intergenic
	};

	
	public static enum Category {
		CDS {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				buf.Item.clear();
				if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
					.sideEffect(r->buf.Item.add(r.toImmutable()))
					.count()>0)
					return true;

				buf.Item.clear();
				if (extraCds!=null && extraCds.ei(hit)
					.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
					.sideEffect(r->buf.Item.add(r.toImmutable()))
					.count()>0)
					return true;
				
				return false;
			}
		},
		RNASEQ,UTR5 {

			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
					
				if (EI.wrap(trs)
							.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
							.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
							.map(tr->tr.getData().get5Utr(tr))
							.filter(tr->tr.getRegion().intersects(hit.getRegion()) && tr.getRegion().isIntronConsistent(hit.getRegion()))
							.sideEffect(r->buf.Item.add(r.toImmutable()))
							.count()>0)
							return true;
					return false;
				}
			
		},
		OffFrame {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				
				if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3!=0)
					.sideEffect(r->buf.Item.add(r.toImmutable()))
					.count()>0)
					return true;

				buf.Item.clear();
				if (extraCds!=null && extraCds.ei(hit)
					.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3!=0)
					.sideEffect(r->buf.Item.add(r.toImmutable()))
					.count()>0)
					return true;
				
				return false;
			}
			
		},
		CDSintoIntron {

			@Override
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (intoIntron!=null && EI.wrap(intoIntron.ei(hit))
						.filter(tr->tr.getData().isCoding())
						.map(tr->tr.getData().getCds(tr))
						.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.induce(hit.getRegion()).getStart()%3==0)
						.sideEffect(r->buf.Item.add(new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), r.getRegion())))
						.count()>0)
						return true;
				return false;
			}
			
		},
		OtherIntoIntron {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (CDSintoIntron.matches(g, hit, extraCds, trs, buf)) return false;
				
				buf.Item.clear();
				if (intoIntron!=null && EI.wrap(intoIntron.ei(hit))
						.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
						.sideEffect(r->buf.Item.add(new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), r.getRegion())))
						.count()>0)
						return true;

				return false;
			}
			
		},Frameshift,
		UTR3 {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (EI.wrap(trs)
						.filter(tr->tr.getData().isCoding() && SequenceUtils.checkCompleteCodingTranscript(g, tr))
						.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
						.map(tr->tr.getData().get3Utr(tr))
						.filter(tr->tr.getRegion().intersects(hit.getRegion()) && tr.getRegion().isIntronConsistent(hit.getRegion()))
						.sideEffect(r->buf.Item.add(r.toImmutable()))
						.count()>0)
						return true;
				return false;
			}
			
		},
		ncRNA {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (EI.wrap(trs)
						.filter(tr->!tr.getData().isCoding() || !SequenceUtils.checkCompleteCodingTranscript(g, tr))
						.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
						.sideEffect(r->buf.Item.add(r.toImmutable()))
						.count()>0)
						return true;
				return false;
			}
			
		},Extra,READS,Substitution,
		Intronic {
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (g.getGenes().ei(hit)
						.filter(gene->gene.contains(hit))
						.count()>0) {
					
					g.getGenes().ei(hit)
						.unfold(gene->g.getTranscripts().ei(gene))
						.unfold(t->EI.wrap(t.getRegion().invert().getParts()))
						.map(p->p.asRegion())
						.filter(p->p.intersects(hit.getRegion()))
						.map(p->new ImmutableReferenceGenomicRegion<Transcript>(hit.getReference(), p))
						.toCollection(buf.Item);
					
						return true;
				}
				return false;
			}
			
		},
		Intergenic {

			@Override
			public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit,
					MemoryIntervalTreeStorage<Transcript> extraCds,
					ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs,
					MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
				if (g.getGenes().ei(hit)
						.filter(gene->gene.contains(hit))
						.count()>0) 
					return false;
				
				ArrayList<GenomicRegion> left = g.getTranscripts().getTree(hit.getReference()).getIntervalsLeftNeighbor(hit.getRegion().getStart(), hit.getRegion().getStart(), new ArrayList<>());
				ArrayList<GenomicRegion> right = g.getTranscripts().getTree(hit.getReference()).getIntervalsRightNeighbor(hit.getRegion().getStop(), hit.getRegion().getStop(), new ArrayList<>());
				if (!left.isEmpty() && !right.isEmpty()) {
					buf.Item.add(new ImmutableReferenceGenomicRegion<>(hit.getReference(), new ArrayGenomicRegion(left.get(0).getEnd(),right.get(0).getStart())));
				}
				
				return true;
			}
			
		},
		PeptideSpliced,AllPeptideSpliced,Unknown;
		
		
		public boolean matches(Genomic g, ImmutableReferenceGenomicRegion<String> hit, MemoryIntervalTreeStorage<Transcript> extraCds, ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs, MutableMonad<List<ImmutableReferenceGenomicRegion<Transcript>>> buf) {
			return false;
		}
	}
	
	public static ImmutableReferenceGenomicRegion<Transcript> getIntron(ImmutableReferenceGenomicRegion<Transcript> t, ImmutableReferenceGenomicRegion<String> loc) {
		ArrayGenomicRegion reg = t.getRegion().invert();
		ArrayGenomicRegion preg = loc.getRegion().intersect(reg);
		return new ImmutableReferenceGenomicRegion<>(t.getReference(), reg.getEnclosingPart(preg.map(preg.getTotalLength()/2)).asRegion(),t.getData());
	}
	
	public static class AnnotationEnum {
		
		private Annotation[] values;
		private HashMap<String,Integer> map;
		private HashMap<String,Integer> catmap;
		public Annotation[] values() {
			return values;
		}
		
		public Annotation valueOf(String name) {
			if (!map.containsKey(name))
				throw new RuntimeException("Annotation "+name+" unknown!");
			return values[map.get(name)];
		}
		
		public Annotation valueOfCategory(String name) {
			if (!catmap.containsKey(name))
				throw new RuntimeException("Category "+name+" unknown!");
			return values[catmap.get(name)];
		}
		
		public boolean hasName(String name) {
			return map.containsKey(name);
		}
		
		private void createMaps() {
			map = ArrayUtils.createIndexMap(values,a->a.name());
			catmap = new HashMap<>();
			for (int i=0; i<values.length; i++)
				for (Category c : values[i].categories)
					catmap.put(c.name(),i);
		}
		
		public boolean contains(Category cat) {
			return EI.wrap(values).filter(a->a.contains(cat)).count()>0;
		}
		
		public HashSet<String> getCategories() {
			return EI.wrap(values()).unfold(a->EI.wrap(a.categories)).map(c->c.name()).set();
		}

		public Annotation getAnnotation(HashSet<String> set) {
			HashSet<Category> cats = EI.wrap(set).map(s->Category.valueOf(s)).set();
			for (Annotation a : values)
				if (a.matches(cats))
					return a;
			throw new RuntimeException(EI.wrap(set).concat(",")+" does not match any annotation!");
		}
		
		private void check() {
			if (EI.wrap(values()).unfold(a->EI.wrap(a.categories)).set().size() != EI.wrap(values()).unfold(a->EI.wrap(a.categories)).list().size())
				throw new RuntimeException("There are categories in multiple annotations!");
		}
		
		public static AnnotationEnum fromParam(boolean all, String cat, boolean noloc) throws IOException {
			if (all) return getAll();
			if (cat.equals("default")) return getDefault(noloc);
			if (cat!=null && new File(cat).exists()) return getFile(cat);
			if (cat!=null) return getDirect(StringUtils.split(cat,','));
			return getDefault(noloc);
		}
		
		public static AnnotationEnum getDefault(boolean noloc) {
			AnnotationEnum re = new AnnotationEnum();
			re.values = EI.wrap(defaultCategories).iff(noloc, ei->ei.chain(EI.singleton(Category.Unknown)))
				.map(c->new Annotation(c.name(),new Category[] {c},re))
				.toArray(Annotation.class);
			re.createMaps();
			return re;
		}

		public static AnnotationEnum getAll() {
			AnnotationEnum re = new AnnotationEnum();
			re.values = EI.wrap(Category.values())
				.map(c->new Annotation(c.name(),new Category[] {c},re))
				.toArray(Annotation.class);
			re.createMaps();
			return re;
		}
		
		public static AnnotationEnum getDirect(String...names) {
			AnnotationEnum re = new AnnotationEnum();
			re.values = EI.wrap(names)
				.map(n->Annotation.parse(re, n))
				.toArray(Annotation.class);
			re.createMaps();
			return re;
		}
		
		public static AnnotationEnum getFile(String file) throws IOException {
			AnnotationEnum re = new AnnotationEnum();
			re.values = EI.lines(file,"#")
				.map(c->Annotation.parse(re,c))
				.toArray(Annotation.class);
			re.createMaps();
			re.check();
			return re;
		}

	}
	
	public static class Annotation {
		
		private String name;
		private Category[] categories;
		private AnnotationEnum aenum;
		
		public Annotation(String name, Category[] categories, AnnotationEnum aenum) {
			this.name = name;
			this.categories = categories;
			this.aenum = aenum;
		}
		
		public boolean matches(Category cat) {
			if (categories.length==1) 
				return categories[0].equals(cat);
			throw new RuntimeException("Orthogonal categories: "+toString()+" - "+cat);
		}

		
		public boolean matches(HashSet<Category> cats) {
			HashSet<Category> c = new HashSet<>(cats);
			c.removeAll(Arrays.asList(categories));
			if (c.size()==0) return true;
			if (c.size()==cats.size()) return false;
			throw new RuntimeException("Orthogonal categories: "+toString()+" - "+EI.wrap(cats).concat(","));
		}

		public boolean contains(Category cat) {
			return EI.wrap(categories).filter(c->c.equals(cat)).count()>0;
		}

		public String name() {
			return name;
		}
		
		public int ordinal() {
			return aenum.map.get(name());
		}
		
		public static Annotation parse(AnnotationEnum aenum, String line) {
			line = line.trim();
			if (line.length()==0) throw new RuntimeException("Cannot parse empty line!");
			String name = null;
			
			int c = line.indexOf(':');
			if (c>0) {
				name = line.substring(0, c).trim();
				line = line.substring(c+1).trim();
			}
			
			String[] cat = StringUtils.split(line, ' ');
			if (name==null) name = cat[0];
			Category[] cats = EI.wrap(cat).map(n->Category.valueOf(n)).toArray(Category.class);
			return new Annotation(name,cats,aenum);
		}
		
		@Override
		public String toString() {
			String re = EI.wrap(categories).map(c->c.name()).concat(" ");
			
			if (!categories[0].name().equals(name))
				re = name+": "+re;

			return re;
		}
	}
	
	public static String removeMod(String s) {
		StringBuilder sb = new StringBuilder(s.length());
		
		int p = 0;
		for (int sepIndex=s.indexOf('('); sepIndex>=0; sepIndex = s.indexOf('(',p)) {
			sb.append(s.substring(p,sepIndex));
			p = s.indexOf(')',p)+1;
		}
		if (sb.length()==0) return s;
		sb.append(s.substring(p));
		return sb.toString();
	}
	
	public static ExtendedIterator<String[]> iteratePeaksLines(String file, MutableMonad<HeaderLine> header, Consumer<String> first) throws IOException {
		LineIterator lit = EI.lines(file);
		String h = lit.next();
		if (header!=null)
			header.Item = new HeaderLine(h,',');
		if (first!=null)
			first.accept(h);
		return lit.map(a->StringUtils.split(a, ','));
	}
	
	public static ExtendedIterator<String[][]> iteratePeaksBlocks(String file, MutableMonad<HeaderLine> header, Consumer<String> first, Predicate<String[]> filterHit, AnnotationEnum restrictToBestAnnotation) throws IOException {
		StringIterator lit = EI.lines(file);
		String h = lit.next();

		if (header==null)
			header = new MutableMonad<>();
		
		header.Item = new HeaderLine(h,',');
		
		int peptide = header.Item.get("Peptide");
		int alc = header.Item.get("ALC (%)",-1);
		if (!header.Item.hasField("Fraction")) {
			lit = lit.mapString(s->s+","+StringUtils.splitField(s, ',', peptide));
			h=h+",Fraction";
			header.Item = new HeaderLine(h,',');
		}
		int fraction = header.Item.get("Fraction");
		int feature = header.Item.hasField("Feature")?header.Item.get("Feature"):-1;
		if (!header.Item.hasField("Scan")) {
			lit = lit.mapString(s->s+","+StringUtils.splitField(s, ',', peptide));
			h=h+",Scan";
			header.Item = new HeaderLine(h,',');
		}
		int scan = header.Item.get("Scan");
		int mode = header.Item.get("mode",-1);
		int id;
		boolean reId = !header.Item.hasField("ID");
		if (!reId)
			id = header.Item.get("ID");
		else {
			id=header.Item.size();
			h = h+",ID";
			header.Item = new HeaderLine(h,',');
		}
		
		if (first!=null)
			first.accept(h);
		
		BiPredicate<String[],String[]> eqFields = (a,b)->{
			return a[id].equals(b[id]);
			
//			if (!a[fraction].equals(b[fraction]))
//				return false;
//			if (!a[scan].equals(b[scan]))
//				return false;
//			if (feature!=-1 && !a[feature].equals(b[feature]))
//				return false;
//			
//			return true;
		};
		
		BitVector bv = new BitVector(64);
		HashSet<String> peps = new HashSet<>();
		
		HashMap<String,String[]> multiCache = new HashMap<>();
		MutablePair<String,String[]> featureCache = new MutablePair<>();
		
		ExtendedIterator<String[][]> re = lit.map(a->splitExtra(a, ',',reId))
				.iff(filterHit!=null, ei->ei.filter(filterHit))
				.iff(reId, ei->ei.sideEffect(a->{
					a[a.length-1] = (a[scan].startsWith("F"+a[fraction]+":")?a[scan]:"F"+a[fraction]+":"+a[scan])+(feature!=-1 && !a[feature].equals("-")?"-"+a[feature]:"");
				}))
				.multiplexUnsorted(eqFields,String[].class)
				.iff(reId, ei->ei.map(m->{ // add multi if necessary
					for (String[] a : m) {
						a[a.length-1] = (a[scan].startsWith("F"+a[fraction]+":")?a[scan]:"F"+a[fraction]+":"+a[scan])+(feature!=-1 && !a[feature].equals("-")?"-"+a[feature]:"");
						if (mode>=0 && a[mode].contains("/")) {
							
							if (feature==-1||"-".equals(a[feature])) {
								String[] x = EI.split(a[scan], ' ').toArray(String.class);
								if (x[0].indexOf(':')!=-1)
									for (int i=1; i<x.length; i++) 
										x[i] = x[0].substring(0, x[0].indexOf(':'))+":"+x[i];
								for (int i=0; i<x.length; i++) {
									multiCache.put(a[fraction]+","+x[i], a);
								}
							}
							else {
								featureCache.set(a[feature], a);
							}
						}
					}
					
					String[] mm;
					if (feature==-1||"-".equals(m[0][feature])) 
						mm = multiCache.remove(m[0][fraction]+","+m[0][scan]);
					else if (m[0][feature].equals(featureCache.Item1)) {
						mm = featureCache.Item2;
					}
					else
						mm = null;
					
					if (mm!=null && ArrayUtils.linearSearch(m, mm)==-1) {
						mm[id] = m[0][id];
						String[][] rm = new String[m.length+1][];
						rm[0] = mm;
						System.arraycopy(m, 0, rm, 1, m.length);
						m = rm;
					}
					if (alc>=0)
						Arrays.sort(m,(a,b)->Integer.compare(Integer.parseInt(b[alc]), Integer.parseInt(a[alc])));
					
					// remove double sequences
					if (bv.size()<m.length)
						bv.setSize(m.length);
					for (int i=0; i<m.length; i++) 
						bv.putQuick(i,peps.add(m[i][peptide]));
					if (bv.cardinality()<m.length)
						m = ArrayUtils.restrict(m, bv);
					
					bv.clear();
					peps.clear();
					return m;
				}));
				
		
		MutableMonad<HeaderLine> uheader = header;
		if (restrictToBestAnnotation!=null) {
			// if there are non-decoys, remove all decoys first!
			re = re.map(af->{
				boolean hasnondecoy = false;
				for (int i=0; i<af.length; i++)
					hasnondecoy |= !af[i][uheader.Item.get("Decoy")].equals("D");
				if (!hasnondecoy)
					return af;
				return ArrayUtils.restrict(af, i->!af[i][uheader.Item.get("Decoy")].equals("D"));
			});
			// now restrict to best category
			re = re.map(af->{
				int best = restrictToBestAnnotation.valueOfCategory(af[0][uheader.Item.get("Category")]).ordinal();
				for (int i=1; i<af.length; i++)
					best = Math.min(best,restrictToBestAnnotation.valueOfCategory(af[i][uheader.Item.get("Category")]).ordinal());

				int ubest = best;
				return ArrayUtils.restrict(af, i->restrictToBestAnnotation.valueOfCategory(af[i][uheader.Item.get("Category")]).ordinal()==ubest);
			});
		}
		
		return re;
	}
	
	private static String[] splitExtra(String s, char separator, boolean extra) {
		if (s.length()==0)
			return new String[0];
		
		ArrayList<String> re = new ArrayList<String>();
		int p = 0;
		for (int sepIndex=s.indexOf(separator); sepIndex>=0; sepIndex = s.substring(p).indexOf(separator)) {
			re.add(s.substring(p,p+sepIndex));
			p += sepIndex+1;
		}
		re.add(s.substring(p));
		if (extra) re.add("");
		return re.toArray(new String[re.size()]);
	}
	
	public static ExtendedIterator<Affinity> iterateNetMHC(String netmhc, String hla, int l, String path, int nthreads) throws IOException {
		// split into 500er packages
		return EI.lines(path).block(500).map(list->{
			try {
				File bfile = File.createTempFile(FileUtils.getNameWithoutExtension(path),".peplist");
				FileUtils.writeAllLines(list, bfile);
				return bfile;
			}catch (IOException e) { throw new RuntimeException(e);}
		}).parallelized(nthreads, 1, ei->ei.unfold(ppath->
			new LineIterator(RunUtils.output(netmhc,"-a",hla,"-inptype","1","-l",l+"","-f",ppath.getPath()))
				//.sideEffect(System.out::println)
				.filter(s->s.length()>0 && !s.startsWith("#") && !s.startsWith("-") && !s.startsWith("Number"))
				.map(s->StringUtils.trim(s))
				.filter(s->!s.toLowerCase().startsWith("error") && !s.toLowerCase().startsWith("pos") && !s.toLowerCase().startsWith("protein") && !s.toLowerCase().contains("distance"))
				.map(s->s.split("\\s+"))
				.map(a->new Affinity(a[1],a[2],Double.parseDouble(a[12])))
				.endAction(()->ppath.delete())
				));
	}

	
	
//	public static class ProcessUnidentifiedProgram extends GediProgram {
//
//
//		public ProcessUnidentifiedProgram(FindGenomicPeptidesParameterSet params) {
//			
//			addInput(params.unidentifiedPeaksOut);
//			addInput(params.hla);
//			addInput(params.netmhc);
//			addInput(params.fdrOut);
//			addInput(params.nthreads);
//			
//			addInput(params.input);
//			
//			// just such that the fdrs are computed
//			addInput(params.pepGenomeFOut);
//			
//			addOutput(params.unidentifiedPeaksFOut);
//		}
//
//		@Override
//		public String execute(GediProgramContext context) throws Exception {
//			
//			File csv = getParameter(0);
//			String hla = getParameter(1);
//			String netmhc = getParameter(2);
//			File fdr = getParameter(3);
//			int nthreads = getIntParameter(4);
//			
//			String prefix = getParameter(5);
//			
//			String input = csv.getPath();
//			
//			
//			HashMap<Integer,double[]> fdrs = new HashMap<>();// len->[score]->[q]
//			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
//			for (String[] a : EI.lines(fdr.getPath())
//				.skip(1, s->header.Item = new HeaderLine(s,','))
//				.map(s->StringUtils.split(s, ','))
//				.filter(s->s[header.Item.get("Mode")].equals("Genome")).loop()) {
//				
//				int len = Integer.parseInt(a[header.Item.get("Peptide length")]);
//				int alc = Integer.parseInt(a[header.Item.get("ALC")]);
//				double q = Double.parseDouble(a[header.Item.get("Q")]);
//				
//				fdrs.computeIfAbsent(len, x->new double[101])[alc] = q;
//			}
//			
//			
//			if (hla==null && new File(prefix+".hla").exists()) 
//				hla = prefix+".hla";
//			
//			if (hla==null && new File(StringUtils.removeFooter(prefix,".csv")+".hla").exists()) 
//				hla = StringUtils.removeFooter(prefix,".csv")+".hla";
//			
//			if (hla!=null) 
//				hla = EI.lines(hla, "#").filter(s->s.length()>0).concat(",");
//			
//			
//			HashMap<String,Affinity> aff = null;
//			
//			if (hla!=null) {
//				try {
//					if (new ProcessBuilder(netmhc,"-h").start().waitFor()==0) {
//						context.getLog().info("netMHC works, will annotate binding affinities! (unidentified)");
//						aff=new HashMap<>();
//					}
//					else
//						context.getLog().info("netMHC does not work, will not annotate binding affinities! (unidentified)");
//				} catch (IOException e) {
//					context.getLog().log(Level.WARNING,"netMHC does not work, will not annotate binding affinities! (unidentified)",e);
//				}
//			} else
//				context.getLog().info("No HLA alleles given, will not annotate binding affinities! (unidentified)");
//			
//			
//			
//			
//			HashMap<String,int[]> best = new HashMap<>();
//			
//			HashMap<Integer,LineOrientedFile> pepout = new HashMap<>();
//			for (String[] a : iteratePeaksLines(input, header,null).loop()) {
//				String pep = removeMod(a[header.Item.get("Peptide")]);
//				
//				int alc = Integer.parseInt(a[header.Item.get("ALC (%)")]);
//				int[] pp = best.computeIfAbsent(pep, x->new int[2]);
//				pp[0] = Math.max(pp[0],alc);
//				pp[1]++;
//				
//				if (aff!=null && pep.length()<=15 && pep.length()>=8) // outside of this range, netMHC does not work!
//					pepout.computeIfAbsent(pep.length(),x->{
//						try {
//							LineOrientedFile re = new LineOrientedFile(File.createTempFile("FindGenomicPeptides.unidentified."+pep.length(), ".peplist").getPath());
//							re.startWriting();
//							return re;
//						} catch (IOException e) {
//							throw new RuntimeException("Could not write to temporary file!",e);
//						}
//					}).writeLine(pep+"\n"+StringUtils.reverse(pep)); // also compute decoys!
//			}
//			for (LineOrientedFile lo : pepout.values()) lo.finishWriting();
//		
//			if (aff!=null)
//				for (Integer l : pepout.keySet()) {
//					context.getLog().info("netMHC predictions for length "+l+" (unidentified)");
//					
//					for (Affinity a : iterateNetMHC(netmhc,hla,l,pepout.get(l).getPath(),nthreads).loop()) {
//						try{
//							aff.merge(a.sequence, a, (x,y)->x.perrank<y.perrank?x:y);
//						} catch (NumberFormatException e) {
//							context.getLog().warning("Unexpected line: "+StringUtils.toString(a));
//						}
//					}
//				}
//			for (LineOrientedFile lo : pepout.values()) lo.delete();
//			
//			
//			
//			context.getLog().info("Annotate list for unidentified");
//
//			HashMap<String, Affinity> uaff = aff;
//			LineWriter out = getOutputWriter(0);
//			for (String[] a : iteratePeaksLines(input, header,
//													h->out.writeLine2(h+",Q,spectra"+(uaff!=null?",HLA allele,netMHC % rank,Decoy HLA allele,Decoy netMHC % rank":"")) 
//													).loop()) {
//				
//				String pep = removeMod(a[header.Item.get("Peptide")]);
//				
//				String affout = "";
//				if (aff!=null) {
//					Affinity af = aff.get(pep);
//					if (af==null)
//						affout=",-,100";
//					else
//						affout=String.format(",%s,%.2f", af.allele,af.perrank);
//					
//					af = aff.get(StringUtils.reverse(pep).toString());
//					if (af==null)
//						affout+=",-,100";
//					else
//						affout+=String.format(",%s,%.2f", af.allele,af.perrank);
//				}
//				
//				
//				boolean thereIsaBetterSpectrum = false;
//				int[] bestHit = best.get(pep);
//				if (bestHit==null || bestHit[0]!=Integer.parseInt(a[header.Item.get("ALC (%)")]))
//					thereIsaBetterSpectrum = true;
//				
//				// remove it from the map (important if there is another spectrum with equal score)!
//				if (!thereIsaBetterSpectrum) 
//					best.remove(pep);
//				
//				
//				if (!thereIsaBetterSpectrum) {
//					
//					int alc = Integer.parseInt(a[header.Item.get("ALC (%)")]);
//					if (fdrs.containsKey(pep.length())) {
//						double q = fdrs.get(pep.length())[alc];
//						
//						out.writef("%s,%.3g,%d%s\n",StringUtils.concat(",", a),q,bestHit[1],affout);
//					}
//				}
//				
//				
//			}
//			
//			
//			out.close();
//
//			try {
//				context.getLog().info("Running R script for plotting");
//				RRunner r = new RRunner(prefix+".mhc.unidentified.R");
//				r.set("file",getOutputFile(0).getPath());
//				r.addSource(getClass().getResourceAsStream("/resources/netMHC.R"));
//				r.run(true);
//			} catch (Throwable e) {
//				context.getLog().log(Level.SEVERE, "Could not plot!", e);
//			}
//			return null;
//		}
//		
//	}
	
	public static class SearchReadsProgram extends GediProgram {

		public SearchReadsProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.pepGenomeAOut);
			addInput(params.reads);
			addInput(params.readsStrandness);
			addInput(params.minlen);
			addInput(params.maxlen);
			addInput(params.nthreads);

			setRunFlag(()->!getParameters(1).isEmpty());
			
			addOutput(params.readSeqTab);
			addOutput(params.readSeqs);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			File csv = getParameter(0);
			String[] readFiles = getParameters(1).isEmpty()?null:getParameters(1).toArray(new String[0]);
			Strandness strandness = getParameter(2);
			int minlen = getIntParameter(3);
			int maxlen = getIntParameter(4);
			int nthreads = getIntParameter(5);

			if (readFiles!=null) {
				if (strandness==null || strandness.equals(Strandness.AutoDetect)) throw new RuntimeException("Specify the strandness of the RNA-seq library!");
			}

			HeaderLine hl = new HeaderLine();
			FilteredIterator<String[]> it = EI.lines(csv).header(hl,',').split(',').filter(a->a[hl.get("Category")].equals(Category.READS.toString()));
			
			context.getLog().info("Collecting peptides identified in category READS ...");
			Trie<SeqList> aho = createTrie(context,hl,it,minlen,maxlen,(a)->new SeqList(a[hl.get("ID")]+","+a[hl.get("Peptide")]));
			context.getLog().info("Preparing for Aho-Corasick");
			aho.prepareAhoCorasick(context.getProgress());
			context.getLog().info("Trie has "+aho.getNodeCount()+" nodes!");

			
			boolean firstForward = !strandness.equals(Strandness.Antisense);
			boolean firstReverse = !strandness.equals(Strandness.Sense);
			
			context.getLog().info("Searching through reads...");
			ExtendedIterator<String[][]> sit = getIter(readFiles);
			
			
			sit.progress(context.getProgress(), -1, r->"Processing reads").parallelized(nthreads, 1024, ()->new HashSet<SeqList>(),(ei,sls)->ei.map(reads->{
				sls.clear();
				if (firstForward && reads[0]!=null) searchAndStoreReads(reads[0][1],aho,sls);
				if (firstReverse && reads[0]!=null) searchAndStoreReads(SequenceUtils.getDnaReverseComplement(reads[0][1]),aho,sls);
				if (firstForward && reads[1]!=null) searchAndStoreReads(reads[1][1],aho,sls);
				if (firstReverse && reads[1]!=null) searchAndStoreReads(SequenceUtils.getDnaReverseComplement(reads[1][1]),aho,sls);
				for (SeqList sl : sls) {
					sl.add(reads);
				}
				return null;
			})).drain();
			
			context.getLog().info("Writing read DNA sequence table ...");
			LineWriter out = getOutputWriter(0);
			out.writeLine("ID,Peptide,DNA Sequence,Read Count");
			for (SeqList sl : aho.values()) {
				for (String dna : sl.pepseqs.elements()) {
					out.writef("%s,%s,%d\n", sl.label,dna,sl.pepseqs.get(dna, 0));
				}
			}
			out.close();
			
			context.getLog().info("Writing read fastqs into "+getOutputFile(1)+"...");
			getOutputFile(1).mkdir();
			for (SeqList sl : aho.values()) {
				LineWriter[] wr = new LineWriter[2];
				if (readFiles.length==2) {
					wr[0] = new LineOrientedFile(getOutputFile(1),sl.label.replace(',', '-').replace(':', '-')+"_R1.fastq.gz").write();
					wr[1] = new LineOrientedFile(getOutputFile(1),sl.label.replace(',', '-').replace(':', '-')+"_R2.fastq.gz").write();
				} else {
					wr[0] = wr[1] = new LineOrientedFile(getOutputFile(1),sl.label.replace(',', '-').replace(':', '-')+".fastq.gz").write();
				}
				
				for (String[][] blocks : sl) {
					if (blocks[0]!=null) EI.wrap(blocks[0]).print(wr[0]);
					if (blocks[1]!=null) EI.wrap(blocks[1]).print(wr[1]);
				}
				
				wr[0].close();
				wr[1].close();
			}
			
			return null;
		}
		@SuppressWarnings("unchecked")
		private ExtendedIterator<String[][]> getIter(String[] readFiles) throws IOException {
			boolean fastq = EI.lines(readFiles[0]).first().startsWith("@");
			if (readFiles.length==1) {
				return EI.lines(readFiles[0])
						.block(fastq?4:2,String.class)
						.map(b->b[0].endsWith("/2")?new String[][] {null,b}:new String[][] {b,null});
			} else if (readFiles.length==2) {
				return EI.lines(readFiles[0])
						.block(fastq?4:2,String.class)
						.fuse(String[].class, EI.lines(readFiles[1])
								.block(fastq?4:2,String.class));
			} else throw new RuntimeException();
		}

		private void searchAndStoreReads(String seq, Trie<SeqList> aho, HashSet<SeqList> re) {
			CharDag dag = new CharDag(CharIterator.fromCharSequence(seq),EI.empty());

			dag.traverse(i->i<3?new TranslateCharDagVisitor<>(aho.ahoCorasickVisitor(l->l*3),SequenceUtils.codeTrie,false):null,(res, varit)->{
				re.add(res.getValue());
				res.getValue().pepseqs.add(seq.substring(res.getStart(),res.getStop()+1));
			});
		}
		public static class SeqList extends Vector<String[][]> {
			private String label;
			private Counter<String> pepseqs = new Counter<String>();

			public SeqList(String label) {
				super();
				this.label = label;
			}
			
		}
	}
	
	
	public static class PlotProgram extends GediProgram {

		public PlotProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.input);
			addInput(params.pepGenomeFOut);
			
			addOutput(params.fdrPlot);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String prefix = getParameter(0);
			
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".fdrstat.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/fdr_stats.R"));
				r.run(true);

			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
			
			return null;
		}
	}
	
	public static class BedOutProgram extends GediProgram {

		public BedOutProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.pepGenomeAOut);
			addInput(params.pepBedOutThreshold);
			addInput(params.pepBedOutMHCThreshold);
			
			addOutput(params.pepBedOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			File csv = getParameter(0);
			double maxq = getDoubleParameter(1);
			double maxmhc = getDoubleParameter(2);
			
			context.getLog().info("Creating bed...");
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 

			LineWriter out = getOutputWriter(0);
			for (String[] a : iteratePeaksLines(csv.getPath(), header,null).loop()) {
				
				double q = header.Item.hasField("Q")?Double.parseDouble(a[header.Item.get("Q")]):0;
				double mhc = 0;
				try {
					if (header.Item.hasField("netMHC % rank"))
						mhc=Double.parseDouble(a[header.Item.get("netMHC % rank")]);
				} catch (NumberFormatException e) {}
				
				boolean decoy = a[header.Item.get("Decoy")].equals("D");
				
				if (q<=maxq && mhc<=maxmhc && !decoy && !a[header.Item.get("Annotation")].equals("Unknown")) {
				
					StringBuilder name = new StringBuilder().append(a[header.Item.get("ID")]).append(" (").append(a[header.Item.get("Category")]).append(" ");
					if (a[header.Item.get("CategoryInfo")].length()>0)
						name.append(a[header.Item.get("CategoryInfo")]).append(",");
					name.append("l=").append(a[header.Item.hasField("Top location count (no decoy)")?header.Item.get("Top location count (no decoy)"):header.Item.get("Location count")]).append(")");
					
					double score = Math.min(100, -Math.log10(q)*100); 

					ImmutableReferenceGenomicRegion<ScoreNameAnnotation> rgr = ImmutableReferenceGenomicRegion.parse(a[header.Item.get("Location")], new ScoreNameAnnotation(name.toString(),score));
					out.writeLine(new BedEntry(rgr).toString());
				}
			}
			out.close();
			
			return null;
		}
	}
	
	
	public static class AnnotateProgram extends GediProgram {


		public AnnotateProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.pepGenomeFOut);
			addInput(params.tis);
			addInput(params.genomic);
			addInput(params.mhcOut);
			addInput(params.rtOut);
			addInput(params.all);
			addInput(params.cat);
			addInput(params.noloc);

			addInput(params.input);
			addInput(params.prefix);
			
			addOutput(params.pepGenomeAOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			String tis = getParameter(1);
			Genomic genomic = getParameter(2);
			File netmhc = getParameter(3);
			File srr = getParameter(4);
			
			boolean call = getBooleanParameter(5);
			String ccat = getParameter(6);
			boolean noloc = getBooleanParameter(7);

			
			String prefix = getParameter(8);
			
			String input = csv.getPath();
			
			
			AnnotationEnum aenum = AnnotationEnum.fromParam(call, ccat, noloc);
			
			context.getLog().info("Annotating: "+EI.wrap(aenum.values()).concat(","));
			if (EI.wrap(aenum.values).filter(ae->ae.matches(Category.CDSintoIntron)).count()>0 ||
					EI.wrap(aenum.values).filter(ae->ae.matches(Category.OtherIntoIntron)).count()>0) {
				if (intoIntron==null) {
					context.getLog().info("Compiling 'into intron' data...");
					int[] intronlen = discoverIntoIntron(genomic,context.getProgress());
					context.getLog().info("Found "+intoIntron.size()+" intron retentions with total intronic lengths of "+intronlen[0]+"nt (CDS) and "+intronlen[1]+"nt (other)!");
				}
			}

			
			GenomicRegionStorage<NameAnnotation> vars = null;
			MemoryIntervalTreeStorage<Transcript> extraCDS = null;
			if (new File(prefix+".var.cit").exists()) {
				context.getLog().info("Reading variants...");
				vars=new CenteredDiskIntervalTreeStorage<NameAnnotation>(new File(prefix+".var.cit").getPath()).toMemory();
				extraCDS = buildFrameshiftCds(genomic,vars);

			}
			GenomicRegionStorage<NameAnnotation> uvars = vars;
			
			Trie<Integer> starts = new Trie<>();
			for (String t : EI.split(tis, ',').map(a->a.toUpperCase().replace('U', 'T')).loop())
				starts.put(t, starts.size());
			
			
			context.getLog().info("Reading netMHC results...");
			HashMap<String,ArrayList<Affinity>> aff = new HashMap<>();
			for (String[] f : EI.lines(netmhc.getPath()).split(',').skip(1).loop()) {
				Affinity a = new Affinity(f[1], f[0], Double.parseDouble(f[2]));
//				aff.merge(a.sequence, a, (x,y)->x.perrank<y.perrank?x:y);
				aff.computeIfAbsent(a.sequence, x->new ArrayList<>()).add(a);
			}
			
			context.getLog().info("Reading SRRCalc results...");
			HashMap<String,String> hi = new HashMap<>();
			for (String[] f : EI.lines(srr).split(',').skip(1).loop()) 
				hi.put(f[0],f[1]);
			
			String intFile = (new File(prefix).getParent()==null?"i":new File(prefix).getParent()+"/i")+new File(prefix).getName();
			
			HashMap<String,String> inti = new HashMap<>();
			if (new File(intFile).exists()) {
				context.getLog().info("Reading intensities...");
				
				MutableMonad<HeaderLine> header = new MutableMonad<>(); 
				for (String[][] f : iteratePeaksBlocks(intFile, header, null,null,null).loop()) {
					inti.put(f[0][header.Item.get("ID")], f[0][header.Item.get("Area")]);
				}
			}
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			RandomNumbers rnd = new RandomNumbers(42);
			
			context.getLog().info("Determining ORFs...");
			
			ArrayList<MutablePair<Integer, List<ImmutableReferenceGenomicRegion<PeptideOrf>>>> selected = new ArrayList<MutablePair<Integer,List<ImmutableReferenceGenomicRegion<PeptideOrf>>>>();
			HashMap<ImmutableReferenceGenomicRegion<Void>,HashSet<String>> orfToPeps = new HashMap<ImmutableReferenceGenomicRegion<Void>, HashSet<String>>();
			
			for (String[] a : iteratePeaksLines(input, header,h->{}).loop()) {
				if (genomic!=null) {
					String[] locs = StringUtils.split(a[header.Item.get("Location")],';');
					String[] catInfo = StringUtils.split(a[header.Item.get("CategoryInfo")],';');
					
					HashSet<String> noinfolocs = new HashSet<String>();
					if (catInfo.length>0)
						for (int i=0; i<locs.length; i++)
							if (catInfo[i].length()==0)
								noinfolocs.add(locs[i]);
					
					BitSet skip = new BitSet(locs.length);
					if (header.Item.hasField("Reads")) {
						// resolve multiple locations by more reads (e.g. if there is one location in the same annotation category with I, and another with L, reads can help
						int[] reads = ArrayUtils.parseIntArray(a[header.Item.get("Reads")],';');
						// set all reads of decoys to -1 (they are only considered if there are decoys only)
						for (int i=0; i<reads.length; i++)
							if (locs[i].startsWith("REV"))
								reads[i] = -1;
						// push all with submaximal reads into skip!
						int max = ArrayUtils.max(reads);
						for (int i=0; i<reads.length; i++) {
							if (reads[i]<max)
								skip.set(i);
						}
					} else {
						// set all decoys to skip
						for (int i=0; i<locs.length; i++)
							if (locs[i].startsWith("REV"))
								skip.set(i);
						if (skip.cardinality()==locs.length)
							skip.clear(); // if there are only decoys
					}
					
					IntArrayList best = new IntArrayList();
					List<List<ImmutableReferenceGenomicRegion<PeptideOrf>>> bestOrfs = new ArrayList<>();
					
					if (!a[header.Item.get("Location")].equals("Extra"))
						for (int i=0; i<locs.length; i++) {
							// skip synonymous mutations!
							if (catInfo.length>0 && catInfo[i].length()>0 && noinfolocs.contains(locs[i]))
								continue;
							
							List<ImmutableReferenceGenomicRegion<PeptideOrf>> orfs = locs[i].equals("")?new ArrayList<>():getOrfs(ImmutableReferenceGenomicRegion.parse(locs[i]), genomic, starts,true);
							orfs.sort((x,y)->x.getData().compareTo(y.getData()));
							
							if (orfs.size()>0 && !locs[i].startsWith("REV") && !skip.get(i)) {
								int cmp = -1;
								if (bestOrfs.isEmpty() || (cmp=orfs.get(0).getData().compareTo(bestOrfs.get(0).get(0).getData()))<=0) {
									if (cmp<0) {
										best.clear();
										bestOrfs.clear();
									}
									best.add(i);
									bestOrfs.add(orfs);  // add just one of the potentially multiple best orfs for this location
								}
							}
						}

					String[] cat = StringUtils.split(a[header.Item.get("Category")],';');
					
					int select = 0;
					List<ImmutableReferenceGenomicRegion<PeptideOrf>> porf = null;
					if (best.isEmpty()) {
						if (skip.cardinality()<locs.length/2) {
							do {
								select=rnd.getUnif(0, locs.length);
							} while (skip.get(select));
						}
						else {
							IntArrayList l = new IntArrayList(locs.length);
							for (int i=0; i<locs.length; i++)
								if (!skip.get(i))
									l.add(i);
							select=l.getInt(rnd.getUnif(0, l.size()));
						}
						
					}
					else {
						int ind = best.size()>1?rnd.getUnif(0, best.size()):0;
						select=best.getInt(ind);
						for (ImmutableReferenceGenomicRegion<PeptideOrf> orf : bestOrfs.get(ind)) { 
							ImmutableReferenceGenomicRegion<Void> lorf = new ImmutableReferenceGenomicRegion<Void>(orf.getReference(),orf.getRegion());
							orfToPeps.computeIfAbsent(lorf, x->new HashSet<>()).add(a[header.Item.get("Location")]);
						}
						porf = bestOrfs.get(ind);
					}
					
					selected.add(new MutablePair<>(select,porf));
				}
					
			}
			
			
			
			context.getLog().info("Annotating list...");
			
			String[] alleles = EI.wrap(aff.values()).unfold(l->EI.wrap(l)).map(a->a.allele).sort().unique(true).toArray(String.class);
			HashMap<String, Integer> alleleIndex = EI.wrap(alleles).indexPosition();
			String alleleHeader = EI.wrap(alleles).concat(",");
			
			HashMap<String,ArrayList<Affinity>> uaff = aff;
			LineWriter out = getOutputWriter(0);
			int index=0;
			for (String[] a : iteratePeaksLines(input, header,
													h->out.writeLine2(h
															+",Hydrophobicity"
															+(inti.size()>0?",Intensity":"")
															+(genomic!=null?",Gene,Symbol,ORF location,ORF TIS prio,ORF TIS,ORF length,ORF nterm,ORF cterm,Region nterm,Region length,Region cterm,Intron index,Upstream sequence,Downstream sequence,ORF Pepcount":"")
															+(genomic!=null && uvars!=null?",Upstream Frameshifting variant":"")
															+(uaff.size()>0?",HLA allele,netMHC % rank,"+alleleHeader:"")) 
													).loop()) {
				
				String hydro = "";
				String intens = "";
				String orf = "";
				String affout = "";
				String fsvar = genomic!=null && vars!=null?",":"";
				
				if (inti.size()>0) {
					String ints = inti.get(a[header.Item.get("ID")]);
					if (ints==null) ints = "NA";
					intens=","+ints;
				}
				
				if (genomic!=null) {
					String[] locs = StringUtils.split(a[header.Item.get("Location")],';');
					String[] catInfo = StringUtils.split(a[header.Item.get("CategoryInfo")],';');
					String[] cat = StringUtils.split(a[header.Item.get("Category")],';');
	
					
					MutablePair<Integer, List<ImmutableReferenceGenomicRegion<PeptideOrf>>> sel = selected.get(index++);
					int select = sel.Item1;
					ImmutableReferenceGenomicRegion<PeptideOrf> orfloc = sel.Item2==null?null:sel.Item2.get(0);
					PeptideOrf porf = sel.Item2==null?null:sel.Item2.get(0).getData();
					int npep = sel.Item2==null?0:EI.wrap(sel.Item2).mapToInt(lorf->orfToPeps.get(new ImmutableReferenceGenomicRegion<>(lorf.getReference(), lorf.getRegion())).size()).max();
					
					String[] seqs = StringUtils.split(a[header.Item.get("Sequence")],';');
					String[] modseq = StringUtils.split(a[header.Item.get("ModSequence")],';');
					a[header.Item.get("Sequence")]=select>=seqs.length?"":seqs[select];
					a[header.Item.get("ModSequence")]=select>=modseq.length?"":modseq[select];
					a[header.Item.get("Location")]=locs.length==0?"":locs[select];
					a[header.Item.get("Category")]=cat[select];
					if (header.Item.hasField("Genome")) {
						String[] genoms = StringUtils.split(a[header.Item.get("Genome")],';');
						a[header.Item.get("Genome")]=genoms.length==0?"":genoms[select];
					}
					a[header.Item.get("CategoryInfo")]=catInfo.length==0?"":catInfo[select];
					if (header.Item.hasField("UniqueCategory"))
						a[header.Item.get("UniqueCategory")]=cat[select];
					
					if (header.Item.hasField("Reads")) {
						String[] reads = StringUtils.split(a[header.Item.get("Reads")],';');
						a[header.Item.get("Reads")]=reads.length==0?"":reads[select];
						String[] dreads = StringUtils.split(a[header.Item.get("Decoyreads")],';');
						a[header.Item.get("Decoyreads")]=dreads.length==0?"":dreads[select];
					}
					
					
					if (porf==null) {
						orf=",-,-,-,-,-,-1,-1,-1,-1,-1,-1,-1,,0,";
					} else {
						ImmutableReferenceGenomicRegion<String> loc = ImmutableReferenceGenomicRegion.parse(locs[select],catInfo.length==0?"":catInfo[select]);
//						String gene = genomic.getGenes().ei(loc).filter(t->t.getRegion().contains(loc.getRegion())).sort((x,y)->{
//							int xl = x.getRegion().intersect(loc.getRegion()).getTotalLength();
//							int yl = y.getRegion().intersect(loc.getRegion()).getTotalLength();
//							int re = Integer.compare(yl, xl);
//							if (re==0) re = Integer.compare(y.getRegion().getTotalLength(), x.getRegion().getTotalLength());
//							return re;
//						}).firstOptional().map(r->r.getData()).orElse("");
//						String symbol = genomic.getGeneTable("symbol").apply(gene);
						
						MutablePair<Category, List<ImmutableReferenceGenomicRegion<Transcript>>> anl = annotateLocation(genomic, aenum, loc, extraCDS);
						
						String gene = EI.wrap(anl.Item2).filter(t->t.getRegion().contains(loc.getRegion())).sort((x,y)->{
							int xl = x.getRegion().intersect(loc.getRegion()).getTotalLength();
							int yl = y.getRegion().intersect(loc.getRegion()).getTotalLength();
							int re = Integer.compare(yl, xl);
							if (re==0) re = Integer.compare(y.getRegion().getTotalLength(), x.getRegion().getTotalLength());
							return re;
						}).firstOptional().map(r->r.getData()==null?null:r.getData().getGeneId()).orElse("");
						
						String frameshift = "";
						if (gene.length()==0) {
							gene = genomic.getGenes().ei(loc).filter(t->t.getRegion().contains(loc.getRegion())).sort((x,y)->{
								int xl = x.getRegion().intersect(loc.getRegion()).getTotalLength();
								int yl = y.getRegion().intersect(loc.getRegion()).getTotalLength();
								int re = Integer.compare(yl, xl);
								if (re==0) re = Integer.compare(y.getRegion().getTotalLength(), x.getRegion().getTotalLength());
								return re;
							}).firstOptional().map(r->r.getData()).orElse("");
						} else {
							int at = gene.indexOf('@');
							if (at>=0) {
								frameshift = gene.substring(at+1);
								gene = gene.substring(0,at);
							}
						}
						
						String symbol = genomic.getGeneTable("symbol").apply(gene);
						
						List<PeptideToRegion> cds = EI.wrap(anl.Item2)
							.map(t->new PeptideToRegion(genomic, t, loc))
							.filter(t->t.region!=null && t.region.getRegion().getTotalLength()>0)
							.list();
						Collections.sort(cds);
						
						int cnterm = -1;
						int ccterm = -1;
						int clen = -1;
						int intronIndex = -1;
						String upstream = "";
						String downstream = "";
						if (cds.size()>0) {
							cnterm = cds.get(0).getNterm();
							ccterm = cds.get(0).getCterm();
							clen = cds.get(0).region.getRegion().getTotalLength();
							intronIndex = cds.get(0).getIntronIndex();
							upstream = cds.get(0).getUpstream(12);
							downstream = cds.get(0).getDownstream(12);
							
//							if (vars!=null && (anl.Item1.equals(Category.CDS) || anl.Item1.equals(Category.OffFrame))) {
//								fsvar+=vars.ei(cds.get(0).getRegion()).filter(vv->cds.get(0).getRegion().induce(vv.getRegion()).getStop()<cds.get(0).getRegion().induce(loc.getRegion()).getStart()).map(vv->{
//									String[] fromto = StringUtils.split(vv.getData().getName(), '>');
//									if ((fromto[1].length()-fromto[0].length())%3!=0) 
//										return vv.toLocationString()+"/"+vv.getData().getName();
//									else
//										return null;
//								}).removeNulls().concat(";");
//							}
						}
						fsvar+=frameshift;
						
						
						orf=String.format(",%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d", gene,symbol,orfloc.toLocationString(),
								porf.getTisPrio(),porf.getTis(),porf.getLength(),porf.getNterm(),porf.getCterm(),
								cnterm,clen,ccterm,
								intronIndex,
								upstream,downstream,
								npep);
					}
				}
				
				
				
				
				String seq = a[header.Item.get("Sequence")];
				hydro = String.format(",%s", hi.getOrDefault(seq, "NA"));

				
				if (aff.size()>0) {
					String[] aout = new String[2+alleles.length];
					Arrays.fill(aout, "100");
					aout[0]="-";
					double best = 100;
					if (aff.get(seq)!=null)
						for (Affinity af : aff.get(seq)) {
							aout[2+alleleIndex.get(af.allele)]=String.format("%.2f",af.perrank);
							if (af.perrank<best) {
								best = af.perrank;
								aout[0]=af.allele;
								aout[1]=String.format("%.2f",af.perrank);
							}
						}
					
					affout = ","+EI.wrap(aout).concat(",");
					
				}
				
				out.writef("%s%s%s%s%s%s\n",StringUtils.concat(",", a),hydro,intens,orf,fsvar,affout);
			}
			
			
			out.close();

			try {
				context.getLog().info("Running R script for plotting");
				RRunner r = new RRunner(prefix+".mhc.R");
				r.set("file",getOutputFile(0).getPath());
				r.addSource(getClass().getResourceAsStream("/resources/netMHC.R"));
				r.run(true);
				
				r = new RRunner(prefix+".rt.R");
				r.set("file",getOutputFile(0).getPath());
				r.addSource(getClass().getResourceAsStream("/resources/hydrophobicity.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}

			return null;
		}
		
	}
	
	public static class NetMHCProgram extends GediProgram {


		public NetMHCProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.pepGenomeFOut);
			addInput(params.hla);
			addInput(params.netmhc);
			addInput(params.nthreads);
			
			addInput(params.input);
			addInput(params.prefix);
			
			addOutput(params.mhcOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			String hla = getParameter(1);
			String netmhc = getParameter(2);
			int nthreads = getIntParameter(3);
			
			String prefix = getParameter(4);
			
			String input = csv.getPath();
			
			if (hla==null && new File(prefix+".hla").exists()) 
				hla = prefix+".hla";
			
			if (hla==null && new File(StringUtils.removeFooter(prefix,".csv")+".hla").exists()) 
				hla = StringUtils.removeFooter(prefix,".csv")+".hla";
			
			if (hla==null && new File(StringUtils.removeFooter(prefix,".csv.gz")+".hla").exists()) 
				hla = StringUtils.removeFooter(prefix,".csv.gz")+".hla";
			
			if (hla!=null) 
				hla = EI.lines(hla, "#").filter(s->s.length()>0).concat(",").replace("*", "");
			
			
			
			HashMap<String,ArrayList<Affinity>> aff = null;
			
			if (hla!=null) {
				try {
					if (new ProcessBuilder(netmhc,"-h").start().waitFor()==0) {
						context.getLog().info("netMHC works, will annotate binding affinities!");
					
						// netMHC predictions
						HashMap<Integer,LineOrientedFile> pepout = new HashMap<>();
						MutableMonad<HeaderLine> header = new MutableMonad<>(); 
						for (String pep : iteratePeaksLines(input, header,null)
								.unfold(a->EI.split(a[header.Item.get("Sequence")], ';').unique(false))
								.loop()) {
							if (pep.length()<=15 && pep.length()>=8) // outside of this range, netMHC does not work!
								pepout.computeIfAbsent(pep.length(),x->{
									try {
										LineOrientedFile re = new LineOrientedFile(File.createTempFile("Prism."+pep.length(), ".peplist").getPath());
										re.startWriting();
										return re;
									} catch (IOException e) {
										throw new RuntimeException("Could not write to temporary file!",e);
									}
								}).writeLine(pep);
						}
						for (LineOrientedFile lo : pepout.values()) lo.finishWriting();
					
						aff = new HashMap<>();
						for (Integer l : pepout.keySet()) {
							context.getLog().info("netMHC predictions for length "+l);
							for (Affinity a : iterateNetMHC(netmhc,hla,l,pepout.get(l).getPath(),nthreads).loop()) {
								try{
									aff.computeIfAbsent(a.sequence, x->new ArrayList<>()).add(a);
								} catch (NumberFormatException e) {
									context.getLog().warning("Unexpected line: "+StringUtils.toString(a));
								}
							}

						}
						for (LineOrientedFile lo : pepout.values()) new File(lo.getPath()).delete();
						
						
					}
					else
						context.getLog().info("netMHC does not work, will not annotate binding affinities!");
				} catch (IOException e) {
					context.getLog().log(Level.WARNING,"netMHC does not work, will not annotate binding affinities!");
				}
			} else
				context.getLog().info("No HLA alleles given, will not annotate binding affinities!");
			
			
			
			context.getLog().info("Writing netMHC predictions...");
			LineWriter out = getOutputWriter(0);
			
			out.writeLine("Sequence,HLA allele,netMHC % rank");
			if (aff!=null)
				for (String s : aff.keySet()) {
					ArrayList<Affinity> l = aff.get(s);
					Collections.sort(l);
					for (Affinity a : l) {
						out.writef("%s,%s,%.2f\n", s,a.allele,a.perrank);
					}
				}
			out.close();
			

			return null;
		}
		
	}
	
	
	public static class SrrCalcProgram extends GediProgram {


		public SrrCalcProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.pepGenomeFOut);
			addInput(params.hydro);
			
			addInput(params.input);
			addInput(params.prefix);
			
			
			addOutput(params.rtOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			boolean hydro = getBooleanParameter(1);
			
			
			LineWriter out = getOutputWriter(0);
			out.writeLine("Sequence,HI");
			
			if (!hydro) {
				context.getLog().info("Skipping SRRCalc (set -hydro to compute).");
			} else {
				String input = csv.getPath();
				context.getLog().info("Querying SRRCalc...");
				Progress progress = context.getProgress().init().setDescription("Contacting webserver (chunks of 500 peptides)");
				MutableMonad<HeaderLine> header = new MutableMonad<>();
				for (String[] peps : iteratePeaksLines(input, header,null)
							.unfold(a->EI.split(a[header.Item.get("Sequence")], ';').unique(false))
							.block(500)
							.map(l->l.toArray(new String[l.size()]))
							.loop()) {
					progress.incrementProgress();
					String[] hi = computeChunk(peps);
					for (int i=0; i<peps.length; i++)
						out.writef("%s,%s\n",peps[i],hi[i]);
				}
				progress.finish();
			}
			out.close();

			return null;
		}
	
		private String[] computeChunk(String[] peps) throws IOException {
			
			String boundary=""+System.currentTimeMillis();
			StringBuilder sb = new StringBuilder();
			sb.append("-----------------------------").append(boundary).append("\r\n");
			sb.append("Content-Disposition: form-data; name=\"seqs\"\r\n\r\n");
			for (String p : peps)
			sb.append(p).append("\r\n");
			sb.append("-----------------------------").append(boundary).append("\r\n");
			sb.append("Content-Disposition: form-data; name=\"sver\"\r\n\r\nssrFAM\r\n-----------------------------").append(boundary).append("\r\n");
			sb.append("Content-Disposition: form-data; name=\"labeld\"\r\n\r\nNONE\r\n-----------------------------").append(boundary).append("\r\n");
			sb.append("Content-Disposition: form-data; name=\"alky\"\r\n\r\nIAM\r\n-----------------------------").append(boundary).append("\r\n");
			sb.append("Content-Disposition: form-data; name=\"chrt\"\r\n\r\nNONE\r\n-----------------------------").append(boundary).append("--\r\n\r\n");



			HttpURLConnection conn = (HttpURLConnection)new URL("http://hs2.proteome.ca/ssrcalc-cgi/SSRCalcQ.pl").openConnection();
			conn.setRequestMethod("POST");
			conn.setRequestProperty("Content-Type","multipart/form-data; boundary=---------------------------"+boundary);
			conn.setDoOutput(true);

			OutputStream os = conn.getOutputStream();
			os.write(sb.toString().getBytes());
			os.flush();
			os.close();
			
			String[] re = new String[peps.length];
			int npep = 0;
			try (LineIterator lit = new LineIterator(conn.getInputStream())) {
				while (lit.hasNext()) {
					String l = lit.next();
					if (l.startsWith("<tr class=\"bodyText\">") && l.indexOf("seq=")!=-1) {
						// check pep order
						String pep = l.substring(l.indexOf("seq=")+4);
						pep = pep.substring(0,pep.indexOf('"'));
						if (!pep.equals(peps[npep])) 
							throw new RuntimeException("Order of peptides in SRRCalc output does not match!");
						String hi = lit.next();
						hi = hi.substring(4,hi.length()-5).trim();
						re[npep++] = hi;
					}
				}
			}				
			if (npep<peps.length)
				 throw new RuntimeException("Insufficient peptides in SRRCalc output!");

			return re;
		}
		
	}
	

	public static class ComputeFdrProgram extends GediProgram {


		public ComputeFdrProgram(FindGenomicPeptidesParameterSet params) {
			
			addInput(params.pepGenomeOut);
			addInput(params.input);
			addInput(params.nthreads);
			addInput(params.all);
			addInput(params.cat);
			addInput(params.genomic);
			addInput(params.noloc);
			
			addInput(params.prefix);
			addOutput(params.pepGenomeFOut);
			addOutput(params.fdrOut);
			
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			
			String in = getParameter(1);
			int nthreads = getIntParameter(2);
			boolean all = getBooleanParameter(3);
			String ccat = getParameter(4);
			Genomic g = getParameter(5);
			boolean noloc = getBooleanParameter(6);
			
			String prefix = getParameter(7);
			
			AnnotationEnum aenum = AnnotationEnum.fromParam(all, ccat, noloc);
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			
			
			String input = csv.getPath();
			
			if (!EI.lines(input).first().contains(",ALC (%)")) {
				context.getLog().info("No ALC score, skipping FDR!");
				
				
				context.getLog().info("Writing peptide FDR list");
				MutablePair<Integer,Integer> key = new MutablePair<>();
				LineWriter out = getOutputWriter(0);
				for (String[][] f : iteratePeaksBlocks(input, header,
														h->out.writeLine2(h+",UniqueCategory,UniqueGenome"), 
														l->!l[header.Item.get("Category")].equals("Unknown"),
														aenum
														).loop()) {
					int len = EI.seq(0, f.length).mapToInt(i->removeMod(f[i][header.Item.get("Peptide")]).length()).unique(false).getUniqueResult(true, true);
					int best = aenum.valueOf(f[0][header.Item.get("Annotation")]).ordinal();
					
					HashSet<String> decoy = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoy")]).set();
					
					String d = "T";
					if (decoy.contains("D"))
						d="D";
					if (decoy.contains("D") && decoy.contains("T"))
						d="B";
					
					key.Item1 = len;
					key.Item2 = best;
					
					key.Item2 = -1;
					
					String se = EI.seq(0, f.length).map(i->f[i][header.Item.get("Sequence")]).concat(";");
					String mse = EI.seq(0, f.length).map(i->f[i][header.Item.get("ModSequence")]).concat(";");
					String geno = EI.seq(0, f.length).map(i->f[i][header.Item.get("Genome")]).concat(";");
					String loc = EI.seq(0, f.length).map(i->f[i][header.Item.get("Location")]).concat(";");
					String cats = EI.seq(0, f.length).map(i->f[i][header.Item.get("Category")]).concat(";"); //.getUniqueResult("More than one annotation!", "No annotation!");
					String catInfo = EI.seq(0, f.length).map(i->f[i][header.Item.get("CategoryInfo")]).concat(";"); 
					String mcats = EI.seq(0, f.length).map(i->Category.valueOf(f[i][header.Item.get("Category")])).sort().unique(true).concat(";");
					String mgenome = EI.seq(0, f.length).map(i->f[i][header.Item.get("Genome")]).sort().unique(true).concat(";");
					
					f[0][header.Item.get("Sequence")]=se;
					f[0][header.Item.get("ModSequence")]=mse;
					f[0][header.Item.get("Genome")]=geno;
					f[0][header.Item.get("Decoy")]=d;
					f[0][header.Item.get("Location")]=loc;
					f[0][header.Item.get("Category")]=cats;
					f[0][header.Item.get("CategoryInfo")]=catInfo;
					
					if (header.Item.hasField("Reads")) {
						String reads = EI.seq(0, f.length).map(i->f[i][header.Item.get("Reads")]).concat(";");
						String dreads = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoyreads")]).concat(";");
						f[0][header.Item.get("Reads")]=reads;
						f[0][header.Item.get("Decoyreads")]=dreads;
					}
					
//					for (int i=0; i<f.length; i++)
						out.writef("%s,%s,%s\n",StringUtils.concat(",", f[0]),mcats,mgenome);
				}
				
				out.close();

				
				LineWriter fout = getOutputWriter(1);
				fout.writeLine("Peptide length,Annotation,ALC,targets,decoys,ambiguous,Q");
				fout.close();
				return null;
			}
			
			
			HashSet<String> allowed = aenum.getCategories();
			
			context.getLog().info("Computing FDR cutoffs");
			
			int[] decoyCounter = {0,0}; // decoy only, both
			int[] decoyCatCounter = new int[aenum.values().length];
			
			
			// [len,prio]->ALC,[target,decoy,both]
			HashMap<MutableTriple<Integer,Integer,String>,int[][]> counter = new HashMap<>();

			for (String[][] f : iteratePeaksBlocks(input, header,null, l->!l[header.Item.get("Category")].equals("Unknown"),aenum).loop()) {
				
				
				int len = EI.seq(0, f.length).mapToInt(i->removeMod(f[i][header.Item.get("Peptide")]).length()).unique(false).getUniqueResult(true, true);
				int alc = EI.seq(0, f.length).mapToInt(i->Integer.parseInt(f[i][header.Item.get("ALC (%)")])).unique(false).getUniqueResult(true, true);
				int cat = aenum.valueOf(f[0][header.Item.get("Annotation")]).ordinal();
		
				HashSet<String> decoy = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoy")]).set();
				if (decoy.contains("D")) 
					decoyCounter[decoy.contains("T")?1:0]++;
				if (decoy.contains("D") && !decoy.contains("T"))
					decoyCatCounter[cat]++;
				
				int ind = 0;
				if (decoy.contains("D"))
					ind=1;
				if (decoy.contains("D") && decoy.contains("T"))
					ind=2;
				
				HashSet<String> genos = EI.wrap(f).map(a->a[header.Item.get("Genome")]).set();
				
				for (String geno : genos)
					counter.computeIfAbsent(new MutableTriple<>(len,cat,geno), x->new int[3][101])[ind][alc]++;
				counter.computeIfAbsent(new MutableTriple<>(len,-1,"-"), x->new int[3][101])[ind][alc]++;
				if (aenum.contains(Category.CDS))
					counter.computeIfAbsent(new MutableTriple<>(len,-2,"-"), x->new int[3][101])[ind][alc]++;
			}
			
//			context.getLog().info("Decoys unique: "+decoyCounter[0]);
//			context.getLog().info("Decoy==target: "+decoyCounter[1]);
//			double fac = 1/(decoyCounter[0]/((double)decoyCounter[0]+decoyCounter[1]));
//			context.getLog().info("Decoy factor: "+fac);
			
			LineWriter oout = new LineOrientedFile(in+prefix+".pep.fdrdata.tsv").write();
			oout.writeLine("Length\tAnnotation\tGenome\tALC\tTarget\tDecoy\tBoth");
			for (MutableTriple<Integer,Integer,String> p : counter.keySet()) {
				int[][] a = counter.get(p);
				for (int i=0; i<=100; i++) {
					oout.writef("%d\t%s\t%s\t%d\t%d\t%d\t%d\n", 
								p.Item1,
								p.Item2==-1?"Total":(p.Item2==-2?"Reference":aenum.values()[p.Item2].name),
								p.Item3,
								i,
								a[0][i],a[1][i],a[2][i]);
				}
			}
			oout.close();
			
			
			
			context.getLog().info("Non parametric fit of mixture score distributions...");
			
			RRunner r = new RRunner(prefix+".fdr.R");
			r.set("prefix",in+prefix);
			r.set("nthreads",nthreads+"");
			r.addSource(getClass().getResourceAsStream("/resources/nonparametric_fit.R"));
			r.run(true);
			
			HashMap<MutableTriple<Integer,Integer,String>,double[][]> qvalPeps = new HashMap<>();
			for (String[] a : EI.lines(in+prefix+".pep.fdr.fit.tsv").header(header.Item).split('\t').loop()) {
				MutableTriple<Integer, Integer,String> p = new MutableTriple<>(Integer.parseInt(a[header.Item.get("Length")]),aenum.hasName(a[header.Item.get("Annotation")])?aenum.valueOf(a[header.Item.get("Annotation")]).ordinal():-1,a[header.Item.get("Genome")]);
				double[][] q = qvalPeps.computeIfAbsent(p, x->new double[101][2]);
				q[Integer.parseInt(a[header.Item.get("ALC")])][0] = a[header.Item.get("FDR")].equals("NA")?1:Double.parseDouble(a[header.Item.get("FDR")]);
				q[Integer.parseInt(a[header.Item.get("ALC")])][1] = a[header.Item.get("PEP")].equals("NA")?1:Double.parseDouble(a[header.Item.get("PEP")]);
			}
			
			
			context.getLog().info("Writing FDR statistics");
			LineWriter fout = getOutputWriter(1);
			fout.writeLine("Peptide length,Annotation,Genome,ALC,targets,decoys,ambiguous,Q");
		
				
			for (MutableTriple<Integer,Integer,String> pair : counter.keySet()) {
				int[][] m = counter.get(pair);
				double[][] a = qvalPeps.get(pair);
				for (int i=0; i<m[0].length; i++) 
					if (pair.Item2>=0)
						fout.writef("%d,%s,%s,%d,%d,%d,%d,%3g\n", pair.Item1, aenum.values()[pair.Item2].name(),pair.Item3, i, m[0][i], m[1][i], m[2][i], a==null?1:a[i][0]);
			}
			fout.close();
			
			double[] q11 = {1,1};
			
			context.getLog().info("Writing peptide FDR list");
			MutableTriple<Integer,Integer,String> key = new MutableTriple<>();
			LineWriter out = getOutputWriter(0);
			for (String[][] f : iteratePeaksBlocks(input, header,
													h->out.writeLine2(h+",UniqueCategory,UniqueGenome,PEP,Q,wrongQ"), 
													null,
													aenum
													).loop()) {
				
				int len = EI.seq(0, f.length).mapToInt(i->removeMod(f[i][header.Item.get("Peptide")]).length()).unique(false).getUniqueResult(true, true);
				int alc = EI.seq(0, f.length).mapToInt(i->Integer.parseInt(f[i][header.Item.get("ALC (%)")])).unique(false).getUniqueResult(true, true);
				String aaa = f[0][header.Item.get("Annotation")];
				int best = aaa.equals("Unknown")?0:aenum.valueOf(aaa).ordinal();
				
				
				HashSet<String> decoy = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoy")]).set();
				
				String d = "T";
				if (decoy.contains("D"))
					d="D";
				if (decoy.contains("D") && decoy.contains("T"))
					d="B";
				if (decoy.contains("U"))
					d="U";
				
				HashSet<String> genos = EI.wrap(f).map(a->a[header.Item.get("Genome")]).set();
					
				String[] genosa = genos.toArray(new String[0]);
				
				int gi = EI.wrap(genosa).mapToDouble(geno->{
					key.Item1 = len;
					key.Item2 = best;
					key.Item3 = geno;
					return (qvalPeps.containsKey(key)?qvalPeps.get(key)[alc]:q11)[0];
				}).argmin();
				
				key.Item1 = len;
				key.Item2 = best;
				key.Item3 = genosa[gi];
				double[] qp = qvalPeps.containsKey(key)?qvalPeps.get(key)[alc]:q11;
				
				key.Item2 = -1;
				double[] qw = qvalPeps.containsKey(key)?qvalPeps.get(key)[alc]:q11;
				
				String se = EI.seq(0, f.length).map(i->f[i][header.Item.get("Sequence")]).concat(";");
				String mse = EI.seq(0, f.length).map(i->f[i][header.Item.get("ModSequence")]).concat(";");
				String loc = EI.seq(0, f.length).map(i->f[i][header.Item.get("Location")]).concat(";");
				String geno = EI.seq(0, f.length).map(i->f[i][header.Item.get("Genome")]).concat(";");
				String cats = EI.seq(0, f.length).map(i->f[i][header.Item.get("Category")]).concat(";"); //.getUniqueResult("More than one annotation!", "No annotation!");
				String catInfo = EI.seq(0, f.length).map(i->f[i][header.Item.get("CategoryInfo")]).concat(";"); 
				String mcats = EI.seq(0, f.length).map(i->Category.valueOf(f[i][header.Item.get("Category")])).sort().unique(true).concat(";");
				String mgenome = EI.seq(0, f.length).map(i->f[i][header.Item.get("Genome")]).sort().unique(true).concat(";");
				
				f[0][header.Item.get("Sequence")]=se;
				f[0][header.Item.get("ModSequence")]=mse;
				f[0][header.Item.get("Genome")]=geno;
				f[0][header.Item.get("Decoy")]=d;
				f[0][header.Item.get("Location")]=loc;
				f[0][header.Item.get("Category")]=cats;
				f[0][header.Item.get("CategoryInfo")]=catInfo;
				
				if (header.Item.hasField("Reads")) {
					String reads = EI.seq(0, f.length).map(i->f[i][header.Item.get("Reads")]).concat(";");
					String dreads = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoyreads")]).concat(";");
					f[0][header.Item.get("Reads")]=reads;
					f[0][header.Item.get("Decoyreads")]=dreads;
				}

				
//				for (int i=0; i<f.length; i++)
					out.writef("%s,%s,%s,%.3g,%.3g,%.3g\n",StringUtils.concat(",", f[0]),mcats,mgenome,qp[1],qp[0],qw[0]);
			}
			
			out.close();

			try {
				context.getLog().info("Running R script for plotting");
				r = new RRunner(prefix+".fdr.R");
				r.set("file",getOutputFile(0).getPath());
				r.addSource(getClass().getResourceAsStream("/resources/fdr_types.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		
			return null;
		}

		
	}
	
	public static String getOrigin(Genomic g, String loc) {
		if (loc.startsWith(readDummyRegion.getReference().getName())) return "READS";
		String o = StringUtils.removeHeader(ImmutableReferenceGenomicRegion.parse(loc).getReference().getName(),"REV_");
		Genomic re = g.getOrigin(o);
		if (re==null) return "-";
		return re.getId();
	}
	
	public static class ToPeptideListProgram extends GediProgram {


		public ToPeptideListProgram(FindGenomicPeptidesParameterSet params) {
			
			addInput(params.annotatedPeaksOut);
			addInput(params.deltaNext);
			addInput(params.deltaFirst);
			addInput(params.all);
			addInput(params.cat);
			addInput(params.noloc);
			addInput(params.writeNextFiltered);
			
			addInput(params.prefix);
			addOutput(params.pepGenomeOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			int minDeltaNext = getIntParameter(1);
			int maxDeltaFirst = getIntParameter(2);
			boolean all = getBooleanParameter(3);
			String cat = getParameter(4);
			boolean noloc = getBooleanParameter(5);
			boolean writeDeltaNext = getBooleanParameter(6);
			
			String input = csv.getPath();
			
			AnnotationEnum aenum = AnnotationEnum.fromParam(all, cat, noloc);
			
			HashSet<String> allowed = aenum.getCategories();
			
			context.getLog().info("Determining best PSM");
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			HashMap<String,int[]> best = new HashMap<>();
			
			for (String[][] f : iteratePeaksBlocks(input, header,null, l->allowed.contains(l[header.Item.get("Category")]),aenum).loop()) {
				for (int i=0; i<f.length; i++) {
					String pep = f[i][header.Item.get("ModSequence")];
					int alc = getScore(header.Item,f[i]);//alcc==-1?-1:Integer.parseInt(f[i][alcc]);
					int[] pp = best.computeIfAbsent(pep, x->new int[2]);
					pp[0] = Math.max(pp[0],alc);
					pp[1]++;
				}
			}
			
			LineWriter out2 = writeDeltaNext?new LineOrientedFile(input+".deltaNext.csv.gz").write():null;
			
			
			context.getLog().info("Writing peptide list");
			LineWriter out = getOutputWriter(0);
			
			for (String[][] f : iteratePeaksBlocks(input, header,
													h->{
														out.writeLine2(h+",Delta next,spectra,Annotation,Top location count,Top location count (no decoy)");
														if (out2!=null)
															out2.writeLine2(h+",Passed deltaNext filter,Passed deltaFirst filter,Passed suboptimal hit for sequence filter");
													},
													l->allowed.contains(l[header.Item.get("Category")]),
													aenum
													).loop()) {
				// identify block with same alc
				int e;
				for (e=1; e<f.length && getScore(header.Item,f[e-1])==getScore(header.Item,f[e]); e++);
				//for (e=1; e<f.length && Integer.parseInt(f[e-1][header.Item.get("ALC (%)")])==Integer.parseInt(f[e][header.Item.get("ALC (%)")]); e++);
				int enodecoy = 0;
				for (int i=0; i<e; i++) {
					if (!f[i][header.Item.get("Location")].startsWith("REV"))
						enodecoy++;
				}
				// this is not reasonable: if there is only the I/L difference between two sequences (i.e. both occurs in the target set of sequences), then why not keep both?
				boolean allSequencesSame = true;
				boolean anySequenceHasBetterSpectrum = false;
				int[] bestHit = null;
				for (int i=0; i<e; i++) {
					if (i>0 && !areEqualUptoIsobar(f[i-1][header.Item.get("ModSequence")],f[i][header.Item.get("ModSequence")]))
						allSequencesSame = false;
					bestHit = best.get(f[i][header.Item.get("ModSequence")]);
					if (bestHit==null || bestHit[0]!=getScore(header.Item,f[i])) 
						anySequenceHasBetterSpectrum = true;
				}
				
				// remove it from the map (important if there is another spectrum with equal score)!
				if (allSequencesSame && !anySequenceHasBetterSpectrum) 
					best.remove(f[0][header.Item.get("ModSequence")]);
				
				int deltaNext = e==f.length?100:(getScore(header.Item,f[0]))-getScore(header.Item,f[e]);
				if (!allSequencesSame) deltaNext=0;
				int deltaFirst = Integer.parseInt(f[0][header.Item.get("Delta first")]);
				
				if (deltaNext>=minDeltaNext && deltaFirst<=maxDeltaFirst && !anySequenceHasBetterSpectrum) {
					
					Annotation anno = aenum.getAnnotation(EI.seq(0,e).map(i->f[i][header.Item.get("Category")]).set());
					for (int i=0; i<e; i++) {
						out.write(StringUtils.concat(",", f[i]));
						out.writef(",%d,%d,%s,%d,%d\n",deltaNext,bestHit[1],anno.name(),e,enodecoy);
					}
				} 
				if (out2!=null) {
					out2.write(StringUtils.concat(",", f[0]));
					out2.writef(",%b,%b,%b\n",deltaNext>=minDeltaNext,deltaFirst<=maxDeltaFirst,!anySequenceHasBetterSpectrum);
				}
				
			}
			
			out.close();
			if (out2!=null)
				out2.close();
			
			return null;
		}
	}

	public static class FindGenomicPeptidesProgram extends GediProgram {



		public FindGenomicPeptidesProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.input);
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.minlen);
			addInput(params.maxlen);
			addInput(params.deltaFirst);
			addInput(params.test);
			addInput(params.rnd);
			addInput(params.all);
			addInput(params.cat);
			addInput(params.writeUnidentified);
			addInput(params.noloc);
			addInput(params.extra);
			addInput(params.rnaseq);
			addInput(params.variants);
			addInput(params.dmod);
			addInput(params.reads);
			addInput(params.readsStrandness);
			
			
			addInput(params.prefix);
			addOutput(params.annotatedPeaksOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String input = getParameter(0);
			Genomic genomic = getParameter(1);
			int nthreads = getIntParameter(2);
			int minlen = getIntParameter(3);
			int maxlen = getIntParameter(4);
			int deltaFirst = getIntParameter(5);
			boolean test = getBooleanParameter(6);
			boolean rnd = getBooleanParameter(7);
			boolean call = getBooleanParameter(8);
			String cat = getParameter(9);
			boolean unident = getBooleanParameter(10);
			boolean noloc = getBooleanParameter(11);
			String extraFasta = getParameter(12);
			String rnaseqFasta = getParameter(13);
			String variants = getParameter(14);
			HashSet<ModifiedAminoAcid> dmod = EI.split(getParameter(15),',').map(s->ModifiedAminoAcid.parseSingle(s)).set();
			String[] reads = getParameters(16).isEmpty()?null:getParameters(16).toArray(new String[0]);
			Strandness strandness = getParameter(17);
			
			if (reads!=null) {
				if (strandness==null || strandness.equals(Strandness.AutoDetect)) throw new RuntimeException("Specify the strandness of the RNA-seq library!");
			}
			
			if (variants==null && new File(input+".var.cit").exists()) 
				variants = input+".var.cit";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv")+".var.cit").exists()) 
				variants = StringUtils.removeFooter(input,".csv")+".var.cit";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv.gz")+".var.cit").exists()) 
				variants = StringUtils.removeFooter(input,".csv.gz")+".var.cit";
			if (variants==null && new File(input+".var.bed").exists()) 
				variants = input+".var.bed";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv")+".var.bed").exists()) 
				variants = StringUtils.removeFooter(input,".csv")+".var.bed";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv.gz")+".var.bed").exists()) 
				variants = StringUtils.removeFooter(input,".csv.gz")+".var.bed";
			if (variants==null && new File(input+".vcf").exists()) 
				variants = input+".vcf";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv")+".vcf").exists()) 
				variants = StringUtils.removeFooter(input,".csv")+".vcf";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv.gz")+".vcf").exists()) 
				variants = StringUtils.removeFooter(input,".csv.gz")+".vcf";
			if (variants==null && new File(input+".var").exists()) 
				variants = input+".var";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv")+".var").exists()) 
				variants = StringUtils.removeFooter(input,".csv")+".var";
			if (variants==null && new File(StringUtils.removeFooter(input,".csv.gz")+".var").exists()) 
				variants = StringUtils.removeFooter(input,".csv.gz")+".var";
			
			
			MutableMonad<RuntimeException> ex = new MutableMonad<RuntimeException>();
			MemoryIntervalTreeStorage<NameAnnotation> vars = null;
			MemoryIntervalTreeStorage<Transcript> extraCDS = null;
			if (variants!=null) {
				context.getLog().info("Reading  "+variants);
				if (variants.endsWith("var"))
					vars = readVar(genomic,variants,ex);
				else if (variants.endsWith("var.cit"))
					vars = readCit(genomic,variants,ex);
				else if (variants.endsWith("var.bed"))
					vars = readBed(genomic,variants,ex);
				else if (variants.endsWith("vcf"))
					vars = readVcf(genomic,variants,ex, context.getLog());
				if (ex.Item!=null) throw ex.Item;
				
				if (vars!=null && !variants.endsWith("var.cit")) {
					new File(input+".var.cit").delete();
					new CenteredDiskIntervalTreeStorage<>(input+".var.cit",NameAnnotation.class).fill(vars);
				}
				
				
				context.getLog().info("Read "+vars.size()+" variants!");
				
				extraCDS = buildFrameshiftCds(genomic,vars);
				context.getLog().info("Frameshifted ORFs: "+extraCDS.size()+" in "+extraCDS.ei().map(g->g.getData().getGeneId()).unique(false).count()+" genes!");
				
			} else
				context.getLog().info("No variant file found!");
			
			
			AnnotationEnum aenum = AnnotationEnum.fromParam(call, cat, noloc);
			
			context.getLog().info("Searching for: "+EI.wrap(aenum.values()).concat(","));
			
			if (EI.wrap(aenum.values).filter(ae->ae.matches(Category.CDSintoIntron)).count()>0 ||
					EI.wrap(aenum.values).filter(ae->ae.matches(Category.OtherIntoIntron)).count()>0) {
				context.getLog().info("Compiling 'into intron' data...");
				int[] intronlen = discoverIntoIntron(genomic,context.getProgress());
				context.getLog().info("Found "+intoIntron.size()+" intron retentions with total intronic lengths of "+intronlen[0]+"nt (CDS) and "+intronlen[1]+"nt (other)!");
			}

			
//			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = new Trie<List<ImmutableReferenceGenomicRegion<String>>>();
//			aho.put("GPYAGKLVAI", new ArrayList<>());
//			MemoryIntervalTreeStorage<NameAnnotation> tvars = new MemoryIntervalTreeStorage<>(NameAnnotation.class);
//			tvars.add(ImmutableReferenceGenomicRegion.parse("test:100-101",new NameAnnotation("C>T")));
//			searchDna("test", genomic.getSequence("3:40499334-40499535").toString(), 0, new RandomNumbers(), aho, 8, 26, tvars);
//			context.getLog().info("Tested!");
			
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = createTrie(context,input,minlen,maxlen,reads!=null);
			
			context.getLog().info("Preparing for Aho-Corasick");
			aho.prepareAhoCorasick(context.getProgress());
			context.getLog().info("Trie has "+aho.getNodeCount()+" nodes!");
			
			
			if (test)
				searchTest(context, genomic,  aho, null,nthreads, minlen,maxlen, rnd, 
						aenum.contains(Category.PeptideSpliced), 
						aenum.contains(Category.Frameshift), 
						aenum.contains(Category.Substitution),vars);
			else
				searchAll(context,genomic,aho,null,nthreads,minlen,maxlen, rnd, 
						aenum.contains(Category.PeptideSpliced), 
						aenum.contains(Category.Frameshift), 
						aenum.contains(Category.Substitution),
						extraFasta,rnaseqFasta,vars);
			
			context.getLog().info("Unifying trie...");
			for (List<ImmutableReferenceGenomicRegion<String>> v : aho.values()) {
				HashSet<ImmutableReferenceGenomicRegion<String>> set = new HashSet<>(v);
				v.clear();
				v.addAll(set);
			}
			
			if (reads!=null) {
				context.getLog().info("Searching through reads...");
				searchAllReads(reads, strandness, rnd?new RandomNumbers():null, aho, nthreads);
			}
			
			context.getLog().info("Writing output");
			
			LineWriter out = getOutputWriter(0);
			LineWriter out2 = unident?new LineOrientedFile(input+".unidentified.csv.gz").write():null;
			ModifiedAminoAcid[] buff = new ModifiedAminoAcid[1024];
			for (int i=0; i<buff.length; i++) buff[i] = new ModifiedAminoAcid();
			int[] ind = new int[maxlen];
			char[] c = new char[maxlen];
			
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			for (String[][] f : iteratePeaksBlocks(input, header, h->{
				out.writeLine2(h+",PSM rank,Location count,Delta first,Decoy,Genome,Location,Category,Sequence,CategoryInfo,ModSequence"+(reads!=null?",Reads,Decoyreads":""));
				if (out2!=null)
					out2.writeLine2(h+",Identified");
			},null,null).loop()) {
				// true adds a new id to each line!
				
				int topalc = -1;
				for (int i=0; i<f.length; i++) {
					
					int alc = getScore(header.Item,f[i]);
					int len = ModifiedAminoAcid.parse(f[i][header.Item.get("Peptide")],buff);
					
					int nloc = 0;
					
					if (len>=minlen && len<=maxlen) {
						if (topalc==-1) topalc = alc;
						
						// count the number of locations
						for (int p=0; p<len; p++)
							c[p] = buff[p].aa;
						do {
							for (int p=0; p<len; p++) {
								if (ind[p]==0)
									c[p] = buff[p].aa;
								else {
									c[p] = buff[p].getAlternative(buff,p);
								}
							}
							String kw = new String(c, 0, len);
							nloc+=aho.get(kw).size();
						}
						while (ArrayUtils.increment(ind, len,(p)->buff[p].getAlternative(buff,p)!=0?2:1));
						
						// and now go through it once more and generate output
						for (int p=0; p<len; p++)
							c[p] = buff[p].aa;
						do {
							for (int p=0; p<len; p++) {
								if (ind[p]==0)
									c[p] = buff[p].aa;
								else {
									c[p] = buff[p].getAlternative(buff,p);
								}
							}
								
							String kw = new String(c, 0, len);
							List<ImmutableReferenceGenomicRegion<String>> add = aho.get(kw);
							if (reads!=null) {
								CountList clist = (CountList)add;
//								if (clist.isEmpty()) { // this can happen, if a peptide is only found in reads.
								// documented out to enable searching with READS as first priority
									if (clist.count.get()>0)
										clist.add(new ImmutableReferenceGenomicRegion<>(readDummyRegion.getReference(), readDummyRegion.getRegion(),kw+","));
									else if (clist.countDecoy.get()>0)
										clist.add(new ImmutableReferenceGenomicRegion<>(readDecoyDummyRegion.getReference(), readDecoyDummyRegion.getRegion(),kw+","));
									// this really adds to the list, and does not count (due to == comparison in CountList)
//								}
							} 
							
							for (ImmutableReferenceGenomicRegion<String> l : add) {
								try {
								Category anno = annotateLocation(genomic, aenum, l, extraCDS).Item1;
								if (aenum.contains(anno) ) {
									out.write(StringUtils.concat(",", f[i]));
									out.writef(",%d,%d,%d,%s,%s,%s,%s,%s,%s",i+1,
											nloc,
											topalc-alc,isDecoy(l)?"D":"T",
											getOrigin(genomic, l.toLocationString()),
											l.toLocationString(),
											anno.toString(),
											l.getData(),
											getModSequence(buff,StringUtils.splitField(l.getData(), ',', 0),dmod));
									
									if (reads!=null) {
										CountList clist = (CountList)add;
										out.writef(",%d,%d",
												clist.count.get(),clist.countDecoy.get());
									}
									out.writeLine();
								}
								} catch (Throwable e) {
									e.printStackTrace();
									System.out.println(l);
									System.exit(1);
								}
//								}
							}
						}
						while (ArrayUtils.increment(ind, len,(p)->buff[p].getAlternative(buff,p)!=0?2:1));
					}
					
					
					if (nloc==0 && noloc) {
						out.write(StringUtils.concat(",", f[i]));
//						out.writef(",,,,,,,,\n");
						len = ModifiedAminoAcid.parse(f[i][header.Item.get("Peptide")],buff);
						StringBuilder seq = new StringBuilder();
						for (int s=0; s<len; s++)
							seq.append(buff[s].aa);
						if (seq.length()>=minlen && seq.length()<=maxlen)
							out.writef(",%d,%d,0,U,%s,,Unknown,%s,,%s\n",
									i+1,
									nloc,
									genomic.getOriginList().get(0),
									seq.toString(),
									getModSequence(buff,f[i][header.Item.get("Peptide")],dmod));
					}
					
					
					if (out2!=null) {
						String pep = removeMod(f[0][header.Item.get("Peptide")]);
						if (/*delta>=deltaNext && */pep.length()>=minlen && pep.length()<=maxlen) {
							out2.write(StringUtils.concat(",", f[i]));
							out2.writef(",%d\n",nloc>0?1:0);
						}
					}
				}
				
				
				
				
			}
			out.close();
			if (out2!=null)
				out2.close();

			
			return null;
		}
		

		private String getModSequence(ModifiedAminoAcid[] pep, String seq, HashSet<ModifiedAminoAcid> dmod) {
			
			StringBuilder sb = new StringBuilder();
			for (int i=0; i<seq.length(); i++) {
				if (dmod.contains(pep[i]))
					sb.append(pep[i].toString());
				else
					sb.append(seq.charAt(i));
			}
			return sb.toString();
		}

	}
	private static MemoryIntervalTreeStorage<NameAnnotation> readCit(Genomic genomic, String variants, MutableMonad<RuntimeException> ex) throws IOException {
		MemoryIntervalTreeStorage<NameAnnotation> re = new CenteredDiskIntervalTreeStorage<NameAnnotation>(variants).toMemory();
		for (ImmutableReferenceGenomicRegion<NameAnnotation> rgr : re.ei().loop()) {
			if (ex.Item==null && !rgr.getData().getName().startsWith("-") && !rgr.getData().getName().startsWith(genomic.getSequence(rgr).toString().toUpperCase()+">")) 
				ex.Item = new RuntimeException("Variant does not match given genome: "+rgr);
		}
		return re;
	}
	
	private static MemoryIntervalTreeStorage<NameAnnotation> readBed(Genomic genomic, String variants, MutableMonad<RuntimeException> ex) throws IOException {
		return Bed.iterateEntries(variants, f->new NameAnnotation(f[3])).sideEffect(rgr->{
			if (ex.Item==null && !rgr.getData().getName().startsWith("-") && !rgr.getData().getName().startsWith(genomic.getSequence(rgr).toString().toUpperCase()+">")) 
				ex.Item = new RuntimeException("Variant does not match given genome: "+rgr);
		}).add(new MemoryIntervalTreeStorage<>(NameAnnotation.class));
	}

	public static MemoryIntervalTreeStorage<NameAnnotation> readVcf(Genomic genomic, String variants,
			MutableMonad<RuntimeException> ex, Logger log) throws IOException {
		Pattern p = Pattern.compile("^(.*):(\\d+)-(\\d+) (.+)$");
		return EI.lines(variants,"#").split('\t').map(a->{
			ImmutableReferenceGenomicRegion<NameAnnotation> rgr = new ImmutableReferenceGenomicRegion<>(
						Chromosome.obtain(a[0]), 
						new ArrayGenomicRegion(Integer.parseInt(a[1])-1,Integer.parseInt(a[1])-1+a[3].length()),
						new NameAnnotation(a[3]+">"+a[4])
						);
			if (rgr.getRegion().isEmpty())
				ex.Item = new RuntimeException("Empty variant: "+StringUtils.concat("\t", a));
			if (ex.Item==null && (genomic.getSequence(rgr)==null || !rgr.getData().getName().startsWith("-") && !rgr.getData().getName().startsWith(genomic.getSequence(rgr).toString().toUpperCase()+">"))) {
				synchronized (log) {
					log.warning("Variant does not match given genome: "+StringUtils.concat("\t", a));	
				}
				return null;
			}
			return rgr;
		}).removeNulls().add(new MemoryIntervalTreeStorage<>(NameAnnotation.class));
	}
	
	private static MemoryIntervalTreeStorage<NameAnnotation> readVar(Genomic genomic, String variants, MutableMonad<RuntimeException> ex) throws IOException {
		Pattern p = Pattern.compile("^(.*):(\\d+)-(\\d+) (.+)$");
		return EI.lines(variants).map(l->{
			Matcher m = p.matcher(l);
			if (m.find()) {
				ImmutableReferenceGenomicRegion<NameAnnotation> rgr = new ImmutableReferenceGenomicRegion<>(
							Chromosome.obtain(m.group(1)), 
							new ArrayGenomicRegion(Integer.parseInt(m.group(2))-1,Integer.parseInt(m.group(3))),
							new NameAnnotation(m.group(4))
							);
				if (rgr.getRegion().isEmpty())
					ex.Item = new RuntimeException("Empty variant: "+l);
				if (ex.Item==null && !rgr.getData().getName().startsWith("-") && !rgr.getData().getName().startsWith(genomic.getSequence(rgr).toString().toUpperCase()+">")) 
					ex.Item = new RuntimeException("Variant does not match given genome: "+l);
				return rgr;
			}
			throw new RuntimeException("Could not parse variant: "+l);
		}).add(new MemoryIntervalTreeStorage<>(NameAnnotation.class));
	}
		

	
	
	private static Trie<List<ImmutableReferenceGenomicRegion<String>>> createTrie(GediProgramContext context, String input, int minlen, int maxlen, boolean reads) throws IOException {
		HeaderLine h = new HeaderLine();
		ExtendedIterator<String[]> lit = EI.lines(input)
				.progress(context.getProgress(), -1, s->s)
				.skip(1,s->h.set(s,','))
				.map(l->StringUtils.split(l, ','));
		
		Function<String[],List<ImmutableReferenceGenomicRegion<String>>> listFac = (a)->reads?new CountList():new Vector<ImmutableReferenceGenomicRegion<String>>();
				
		return createTrie(context, h, lit, minlen, maxlen, listFac);
	}
	
	private static <T> Trie<T> createTrie(GediProgramContext context, HeaderLine h, ExtendedIterator<String[]> lit, int minlen, int maxlen, Function<String[],T> listFac) throws IOException {
		context.getLog().info("Creating keyword trie");
		
		ModifiedAminoAcid[] buff = new ModifiedAminoAcid[1024];
		for (int i=0; i<buff.length; i++) buff[i] = new ModifiedAminoAcid();
		int[] ind = new int[maxlen];
		char[] c = new char[maxlen];
		
		Trie<T> aho = new Trie<>();
		for (String[] a : lit.loop()) {
		
			String l = a[h.get("Peptide")];
			int len = ModifiedAminoAcid.parse(l,buff);
			if (len<minlen || len>maxlen) // returned -1 if >buff (i.e. maxlen)
				continue;
			
			
			for (int i=0; i<len; i++)
				c[i] = buff[i].aa;
			do {
				for (int i=0; i<len; i++) {
					if (ind[i]==0)
						c[i] = buff[i].aa;
					else {
						c[i] = buff[i].getAlternative(buff,i);
					}
				}
					
				String kw = new String(c, 0, len);
				if (!aho.containsKey(kw)) {
					aho.put(kw, listFac.apply(a));
				}
				
//					if (l.startsWith("Q(-17.03)")) {
//						c[0] = 'E';
//						String kw = new String(c, 0, len);
//						if (!aho.containsKey(kw))
//							aho.put(kw, list); // this might get overwritten, if this peptide also occurs (this is handled above when annotating the peaks file)
//					}
				
			}
			while (ArrayUtils.increment(ind, len,(i)->buff[i].getAlternative(buff,i)!=0?2:1));
		}
		
		context.getLog().info("...contains "+aho.size()+" words!");
		
		
		
		return aho;
	}
	
	public static class CountList extends Vector<ImmutableReferenceGenomicRegion<String>> {
		private AtomicInteger count = new AtomicInteger();
		private AtomicInteger countDecoy = new AtomicInteger();

		@Override
		public boolean add(ImmutableReferenceGenomicRegion<String> e) {
			if (e==readDecoyDummyRegion) {
				countDecoy.addAndGet(1);
				return true;
			}
			else if (e==readDummyRegion) {
				count.addAndGet(1);
				return true;
			}
			return super.add(e);
		}
		
	}
	
	
	public static class PeptidesCountProgram extends GediProgram {



		public PeptidesCountProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.genomic);
			addInput(params.nthreads);
			addInput(params.test);
			addInput(params.cat);
			addInput(params.extra);
			addInput(params.rnaseq);
			addInput(params.all);
			addInput(params.anchors);
			
			addInput(params.prefix);
			addOutput(params.countOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			Genomic genomic = getParameter(0);
			int nthreads = getIntParameter(1);
			boolean test = getBooleanParameter(2);
			String cat = getParameter(3);
			String extraFasta = getParameter(4);
			String rnaseqFasta = getParameter(5);
			boolean call = getBooleanParameter(6);
			String anchors = getParameter(7);
			
			MemoryIntervalTreeStorage<Transcript> extraCDS = null;
			
			
			AnnotationEnum aenum = AnnotationEnum.fromParam(call, cat, false);
			context.getLog().info("Searching for: "+EI.wrap(aenum.values()).concat(","));
			
			
			context.getLog().info("Preparing peptide counter");
			PeptideCounter counter = new PeptideCounter(context,anchors,aenum.values.length, ahit -> aenum.valueOf(annotateLocation(genomic,aenum, ahit, extraCDS).Item1.name()).ordinal());
			
			context.getLog().info("Counter has "+counter.length()+" entries, "+(long)counter.anno.length()*Long.BYTES+" bytes!");
			
			genomic.getGenes();

			
			if (aenum.contains(Category.AllPeptideSpliced)) {
				String[] allprot = genomic.getTranscripts().ei().filter(r->r.getData().isCoding())
						.map(r->SequenceUtils.translate(genomic.getSequence(r.getData().getCds(r)).toString().toUpperCase()))
						.map(s->StringUtils.removeFooter(s, "*").replace('I', 'L'))
						.iff(test,ei->ei.head(1000))
						.toArray(String.class);
				counter.countSpliced(allprot, aenum.valueOf("AllPeptideSpliced").ordinal(),context,nthreads);
			}
			

			if (test)
				searchTest(context, genomic,  null,counter, nthreads, 9,9, false, 
						aenum.contains(Category.PeptideSpliced), 
						aenum.contains(Category.Frameshift), 
						aenum.contains(Category.Substitution),null);
			else
				searchAll(context,genomic,null,counter,nthreads,9,9, false, 
						aenum.contains(Category.PeptideSpliced), 
						aenum.contains(Category.Frameshift), 
						aenum.contains(Category.Substitution),
						extraFasta,rnaseqFasta,null);
			long[] table = counter.getTable();
			
			LineWriter out = getOutputWriter(0);
			out.writeLine("Anchor\tCategory\tCount");
			for (Annotation a : aenum.values)
				out.writef("%s\t%s\t%d\n",anchors,a.name(),table[a.ordinal()]);
			out.writef("%s\tNone\t%d\n",anchors,table[(int) counter.bitmask]);
			out.close();

			
			return null;
		}

	}
	
	private static class PeptideCounter  {
		static final int LEN = 9;
		
		private ToIntFunction<ImmutableReferenceGenomicRegion<String>> annotator;
		private int[][] alphas = new int[LEN][];
		private int[] alphaSize = new int[LEN];
		private AtomicLongArray anno;
		private int[] anchorPos;
		private int[] nonAnchorPos;
		private long total;
		
		private int bits;
		private long bitmask;
		private int blocks;

		public PeptideCounter(GediProgramContext context, String anchors, int categories,ToIntFunction<ImmutableReferenceGenomicRegion<String>> annotator) {
			this.annotator = annotator;
			
			categories++;
			bits = (int)Math.ceil(Math.log(categories)/Math.log(2));
			bitmask = (1<<bits)-1;
			blocks=Long.BYTES*8/bits;
			
			// 1:LVTM,8:KYR
			// 1:LM,8:VL
			char[] aa19 = EI.wrap(SequenceUtils.code.values()).filter(aa->!aa.equals("I")).filter(aa->!aa.equals("*")).unique(false).concat().toCharArray();
			Arrays.sort(aa19);
			if (aa19.length!=19) throw new RuntimeException(String.valueOf(aa19));
			
			int[] alpha19 = new Alphabet(aa19).createIndex(); 
			Arrays.fill(alphas, alpha19);
			Arrays.fill(alphaSize, aa19.length);
			
			for (String a : EI.split(anchors, ',').loop()) {
				String[] b = StringUtils.split(a, ":");
				if (b.length!=2) throw new RuntimeException("Anchor description has wrong format!");
				int p=Integer.parseInt(b[0]);
				char[] combi = EI.substrings(b[1],1).filter(aa->alpha19[aa.charAt(0)]!=-1).unique(false).concat().toCharArray();
				context.logf("Anchor residue %d: %s", p, String.valueOf(combi));
				this.alphas[p] = new Alphabet(combi).createIndex();
				this.alphaSize[p] = combi.length;
			}
			
			anchorPos = EI.seq(0, this.alphas.length).filterInt(i->this.alphas[i]!=alpha19).toIntArray();
			nonAnchorPos = EI.seq(0, this.alphas.length).filterInt(i->this.alphas[i]==alpha19).toIntArray();
			
			
			total = 1;
			for (int i=0; i<this.alphas.length; i++)
					total*=alphaSize[i];
			if (total/blocks>Integer.MAX_VALUE) throw new RuntimeException(total/blocks+">"+Integer.MAX_VALUE);
			
			long[] a = new long[(int) ((total+blocks-1)/blocks)];
			long allmax = 0;
			for (int i=0; i<blocks; i++)
				allmax|=bitmask<<(i*bits);
			Arrays.fill(a, allmax);
			
			anno = new AtomicLongArray(0);// byte[(int)total];
			try {
				ReflectionUtils.set(anno,"array",a);
			} catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
				e.printStackTrace();
			}

		}

		public long length() {
			return total;
		}

		public long[] getTable() {
			long[] table = new long[(int) (bitmask+1)];
			
			for (long i=0; i<total; i++) {
				int apos = (int) (i/blocks);
				int offset = (int) ((i%blocks)*bits);
				long v = anno.get(apos);
				table[(int) ((v>>>(offset))&bitmask)]++;
			}
			
			return table;
		}

		public void countSpliced(CharSequence[] sequences, int category, GediProgramContext context, int nthreads) {
			for (int cut=1; cut<LEN; cut++) {
				long num = countSpliced(sequences, cut,category,context,nthreads);
//				System.out.println(cut+" "+num);
			}
		}
		public long countSpliced(CharSequence[] sequences, int cut, int category, GediProgramContext context, int nthreads) {
			HashSet<Long> left = new HashSet<Long>();
			for (CharSequence s : sequences) 
				EI.substrings(s.toString(), cut, true).map(part->{
					for (int i=0; i<anchorPos.length && anchorPos[i]<part.length(); i++)
						if (alphas[anchorPos[i]][part.charAt(anchorPos[i])]==-1)
							return null;
					long index = 0;
					for (int i=0; i<part.length(); i++) {
						int charIndex = alphas[i][part.charAt(i)];
						if (charIndex==-1) 
							return null;
						index=index*alphaSize[i]+charIndex;
					}
					for (int i=part.length(); i<LEN; i++)
						index=index*alphaSize[i];
					return index;
				}).removeNulls().toCollection(left);
			
			HashSet<Integer> right = new HashSet<Integer>();
			for (CharSequence s : sequences) 
				EI.substrings(s.toString(), LEN-cut, true).map(part->{
					for (int i=anchorPos.length-1; i>=0 && anchorPos[i]-cut>=0; i--)
						if (alphas[anchorPos[i]][part.charAt(anchorPos[i]-cut)]==-1)
							return null;
					int index = 0;
					for (int i=0; i<part.length(); i++) {
						int charIndex = alphas[i+cut][part.charAt(i)];
						if (charIndex==-1) 
							return null;
						index=index*alphaSize[i+cut]+charIndex;
					}
					return index;
				}).removeNulls().toCollection(right);

//			EI.wrap(sequences).progress(context.getProgress(), sequences.length, s->longerPart+"+"+(LEN-longerPart)).parallelized(nthreads, 10, ei->ei.map(s->{
//				
//				for (String longer : EI.substrings(s.toString(), longerPart,true).loop()) {
//					for (String other : shorter) {
//						add(longer,other,category);
//						add(other,longer,category);
//					}
//				}
//				return true;
//			})).drain();
			
			long[] raw;
			try {
				raw = ReflectionUtils.get(anno, "array");
			} catch (NoSuchFieldException | SecurityException | IllegalArgumentException | IllegalAccessException e) {
				throw new RuntimeException(e);
			}
			
			EI.wrap(left).progress(context.getProgress(), left.size(), x->"Cut point: "+cut).parallelized(nthreads, 1, ei->ei.map(l->{
				for (Integer r : right) {
					add(l.longValue()+r.longValue(),category);
				}
				return null;
			})).drain();
			
//			Progress progress = context.getProgress().init();
//			progress.setDescription("Cut point: "+cut);
//			progress.setCount(left.size());
//			for (Long l : left) {
//				progress.incrementProgress();
//				for (Integer r : right) {
//					add(l.longValue()+r.longValue(),category,raw);
//				}
//			}
//			progress.finish();
			
			return left.size()*(long)right.size();
		}
		
		
		private void add(long index, int val) {
			int apos = (int) (index/blocks);
			int offset = (int) ((index%blocks)*bits);
			anno.accumulateAndGet(apos, val, (ori,n)->{
				long current = (ori>>>offset)&bitmask;
				if (current<=n) return ori;
				return (ori&~(bitmask<<offset)) | n<<offset;
			});
		}
		
		private void add(long index, int val, long[] a) {
			int apos = (int) (index/blocks);
			int offset = (int) ((index%blocks)*bits);
			long ori = a[apos];
			a[apos]=(ori&~(bitmask<<offset)) | val<<offset;
		}

		public boolean count(CharSequence sequence, BiFunction<Integer,String,ImmutableReferenceGenomicRegion<String>> e) {
			final int N=sequence.length()-LEN+1; 
			seqloop:for (int s=0; s<N; s++) {
				for (int i=0; i<anchorPos.length; i++)
					if (alphas[anchorPos[i]][sequence.charAt(s+anchorPos[i])]==-1)
						continue seqloop;
				
				long index = 0;
				for (int i=0; i<LEN; i++) {
					int charIndex = alphas[i][sequence.charAt(s+i)];
					if (charIndex==-1) 
						continue seqloop;
					index=index*alphaSize[i]+charIndex;
				}
				ImmutableReferenceGenomicRegion<String> rgr = e.apply(s,sequence.subSequence(s, s+LEN).toString());
				if (rgr!=null) {
					int apos = (int) (index/blocks);
					int offset = (int) ((index%blocks)*bits);
					int val = annotator.applyAsInt(rgr);
//					decodeAndPrint(anno.get(apos));
//					System.out.println(val+" "+offset/bits);
					anno.accumulateAndGet(apos, val, (ori,n)->{
						long current = (ori>>>offset)&bitmask;
						if (current<=n) return ori;
						return (ori&~(bitmask<<offset)) | n<<offset;
					});
				}
//					anno[index] = (byte) Math.min(annotator.applyAsInt(rgr),anno[index]);
			}
			
			
			return true;
		}

		private void decodeAndPrint(long v) {
			StringBuilder sb = new StringBuilder();
			for (int o=0; o<blocks; o++) 
				sb.append((int) ((v>>(o*bits))&bitmask)).append(",");
			System.out.println(sb.substring(0,sb.length()-1));
		}
		private void decodeAndCheck(long v) {
			for (int o=0; o<blocks; o++) {
				int vv = (int) ((v>>(o*bits))&bitmask);
				if (vv>8 && vv!=bitmask) {
					decodeAndPrint(v);
					return;
				}
			}
		}
		
	}
	
	private static void searchAll(GediProgramContext context, Genomic genomic, 
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, PeptideCounter counter, int nthreads, int minlen, int maxlen, boolean rnd, boolean spliced, boolean frameshift, boolean substitution, String extraFasta, String rnaseqFasta, MemoryIntervalTreeStorage<NameAnnotation> vars) {			
		
		String[] fasta = genomic.getGenomicFastaFiles().toArray(String.class);
		FastaHeaderParser pars = new DefaultFastaHeaderParser(' ');
		HeaderLine mh = new HeaderLine();
		Pattern mp = Pattern.compile("^p\\.([A-Z])\\d+([A-Z])$");
		
		
//			String sss = genomic.getSequence("X:0-"+genomic.getLength("X")).toString();
//			EI.seq(0, 12).map(frame->new MutablePair<>("", frame)).parallelized(nthreads, 1, ()->(rnd?new RandomNumbers():null), (ei,ran)->ei
//					.map(ss->{
//						String seq = sss;
//						String name = "X";
//						int frame = ss.Item2;
//						if (isDna(seq) || (frame>=202 && frame<=207))
//							return searchDna(name,seq,frame,ran,aho,minlen,maxlen,vars);
//						else 
//							return searchProtein(name,seq,frame,ran,aho);
//					}))
//					.iff(genomic==null, ei->ei.map(a->null))
//					.removeNulls()
//					.log(context.getLog(),Level.INFO);
//			
//			for (List<ImmutableReferenceGenomicRegion<String>> v : aho.values()) {
//				HashSet<ImmutableReferenceGenomicRegion<String>> set = new HashSet<>(v);
//				v.clear();
//				v.addAll(set);
//			}
//			
//			if (true) return;
		
		context.getLog().info("Starting Aho-Corasick using "+nthreads+" threads");
		EI.wrap(fasta).unfold(s->new FastaFile(s).entryIterator2())
			.iff(genomic==null,ei->ei.progress(context.getProgress(), -1, fe->"Processing: "+fe.getHeader()))
			.sideEffect(fe->{	
				if (pars.getId(fe.getHeader()).contains(":") && genomic!=null)
					throw new RuntimeException("Sequence names may not contain : !");
			})
			.iff(genomic!=null, ei->
				genomic.getTranscripts().ei()
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Peptides over known exon-exon junctions...")
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.chain(ei))
//				.head(1)
			.unfold(fe->{
				if (isDna(fe.getSequence()))
					return EI.seq(0, 12).map(frame->new MutablePair<>(fe, frame));
				else
					return EI.seq(0, 2).map(frame->new MutablePair<>(fe, frame));
			})
			.iff(extraFasta!=null, ei->
			EI.wrap(new FastaFile(extraFasta)).unfold(ff->ff.entryIterator2())
				.progress(context.getProgress(), -1, x->"Extra fasta...")
				.map(tr->new FastaEntry(pars.getId(tr.getHeader()), tr.getSequence()))
				.filter(fe->fe.getSequence().length()>=minlen)
				.unfold(fe->EI.wrap(new MutablePair<>(fe,102),new MutablePair<>(fe,103)))
				.chain(ei))
			.iff(genomic!=null && spliced, ei->
				genomic.getTranscripts().ei()
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Spliced peptides...")
					.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
					.map(tr->tr.getData().getCds(tr))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.unfold(fe->EI.wrap(new MutablePair<>(fe,42),new MutablePair<>(fe,43)))
					.chain(ei))
			.iff(genomic!=null && frameshift, ei->
				genomic.getTranscripts().ei()
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Frameshift peptides...")
					.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
					.map(tr->tr.getData().getCds(tr).alterRegion(cds->cds.union(tr.getData().get3Utr(tr).getRegion())))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.unfold(fe->EI.wrap(new MutablePair<>(fe,84),new MutablePair<>(fe,85)))
					.chain(ei))
			.iff(genomic!=null && substitution, ei->
				genomic.getTranscripts().ei()
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Substitution peptides...")
					.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
					.map(tr->tr.getData().getCds(tr))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.unfold(fe->EI.wrap(new MutablePair<>(fe,92),new MutablePair<>(fe,93)))
					.chain(ei))
			.iff(rnaseqFasta!=null, ei->
					EI.wrap(new FastaFile(rnaseqFasta)).unfold(ff->ff.entryIterator2())
						.progress(context.getProgress(), -1, x->"RNA-seq fasta...")
						.map(tr->new FastaEntry(pars.getId(tr.getHeader()).replace(':','_'), tr.getSequence()))
						.filter(fe->fe.getSequence().length()>=minlen)
						.unfold(fe->EI.wrap(new MutablePair<>(fe,202),new MutablePair<>(fe,203),new MutablePair<>(fe,204),new MutablePair<>(fe,205),new MutablePair<>(fe,206),new MutablePair<>(fe,207)))
						.chain(ei))
			.parallelized(nthreads, 1, ()->(rnd?new RandomNumbers():null), (ei,ran)->ei
			.map(ss->{
				String seq = ss.Item1.getSequence();
				String name = pars.getId(ss.Item1.getHeader());
				int frame = ss.Item2;
				
				if (isDna(seq) || (frame>=202 && frame<=207))
					return counter!=null?countDna(name, seq, frame, counter, minlen, maxlen): searchDna(name,seq,frame,ran,aho,minlen,maxlen,vars);
				else 
					return counter!=null?countProtein(name, seq, frame, counter):searchProtein(name,seq,frame,ran,aho);

			}))
			.iff(genomic==null, ei->ei.map(a->null))
			.removeNulls()
			.log(context.getLog(),Level.INFO);
		
		
	}
	
	private static void searchTest(GediProgramContext context, Genomic genomic, 
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, PeptideCounter counter, int nthreads, int minlen, int maxlen, boolean rnd, boolean spliced, boolean frameshift, boolean substitution, MemoryIntervalTreeStorage<NameAnnotation> vars) {
		
		String chr = "18";
		
		FastaHeaderParser pars = new DefaultFastaHeaderParser(' ');
		
		String[] fasta = genomic.getGenomicFastaFiles().toArray(String.class);
		
//		EI.singleton(new MutablePair<>(new FastaEntry("9",genomic.getSequence(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain("9+"), new ArrayGenomicRegion(0,(int)genomic.getLength("9")))).toString()),1))
//			.map(ss->{
//				String seq = ss.Item1.getSequence();
//				String name = pars.getId(ss.Item1.getHeader());
//				int frame = ss.Item2;
////				if (isDna(seq))
////					return searchDna(name,seq,frame,new RandomNumbers(),aho,minlen,maxlen,vars);
////				else 
////					return searchProtein(name,seq,frame,new RandomNumbers(),aho);
//				if (isDna(seq))
//					return counter!=null?countDna(name, seq, frame, counter, minlen, maxlen):searchDna(name,seq,frame,new RandomNumbers(),aho,minlen,maxlen,vars);
//				else 
//					return counter!=null?countProtein(name, seq, frame, counter):searchProtein(name,seq,frame,new RandomNumbers(),aho);
//			})
//			.iff(genomic==null, ei->ei.map(a->null))
//			.removeNulls()
//			.log(context.getLog(),Level.INFO);
		
		context.getLog().info("Starting Aho-Corasick using "+nthreads+" threads");
		context.getLog().info("Only chromosome "+chr);
		EI.singleton(new FastaEntry(chr,genomic.getSequence(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(chr+"+"), new ArrayGenomicRegion(0,(int)genomic.getLength(chr)))).toString()))
			.sideEffect(fe->{
				if (pars.getId(fe.getHeader()).contains(":") && genomic!=null)
					throw new RuntimeException("Sequence names may not contain : !");
			})
			.iff(genomic!=null, ei->
				genomic.getTranscripts().ei(chr+"+").chain(genomic.getTranscripts().ei(chr+"-"))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.chain(ei))
			.unfold(fe->{
				if (isDna(fe.getSequence()))
					return EI.seq(0, 12).map(frame->new MutablePair<>(fe, frame));
				else
					return EI.seq(0, 1).map(frame->new MutablePair<>(fe, frame));
			})
			.iff(genomic!=null && spliced, ei->
				genomic.getTranscripts().ei(chr+"+")
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Spliced peptides...")
					.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
					.map(tr->tr.getData().getCds(tr))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.unfold(fe->EI.wrap(new MutablePair<>(fe,42),new MutablePair<>(fe,43)))
					.chain(ei))
			.iff(genomic!=null && frameshift, ei->
			genomic.getTranscripts().ei(chr+"+")
				.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Frameshift peptides...")
				.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
				.map(tr->tr.getData().getCds(tr).alterRegion(cds->cds.union(tr.getData().get3Utr(tr).getRegion())))
				.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
				.filter(fe->fe.getSequence().length()>=minlen*3)
				.unfold(fe->EI.wrap(new MutablePair<>(fe,84),new MutablePair<>(fe,85)))
				.chain(ei))
			.iff(genomic!=null && substitution, ei->
				genomic.getTranscripts().ei(chr+"+")
					.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Substitution peptides...")
					.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
					.map(tr->tr.getData().getCds(tr))
					.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
					.filter(fe->fe.getSequence().length()>=minlen*3)
					.unfold(fe->EI.wrap(new MutablePair<>(fe,92),new MutablePair<>(fe,93)))
					.chain(ei))
			.parallelized(nthreads, 1, ()->(rnd?new RandomNumbers():null), (ei,ran)->ei
			.map(ss->{
				String seq = ss.Item1.getSequence();
				String name = pars.getId(ss.Item1.getHeader());
				int frame = ss.Item2;
				if (isDna(seq))
					return counter!=null?countDna(name, seq, frame, counter, minlen, maxlen):searchDna(name,seq,frame,ran,aho,minlen,maxlen,vars);
				else 
					return counter!=null?countProtein(name, seq, frame, counter):searchProtein(name,seq,frame,ran,aho);
			}))
			.iff(genomic==null, ei->ei.map(a->null))
			.removeNulls()
			.log(context.getLog(),Level.INFO);
		
	}

	private static String searchAllReads(String[] files, Strandness strandness, RandomNumbers rnd,
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, int nthreads) throws IOException {

		boolean firstForward = !strandness.equals(Strandness.Antisense);
		boolean firstReverse = !strandness.equals(Strandness.Sense);
//		boolean secondForward = files.length==2 && firstReverse;
//		boolean secondReverse = files.length==2 && firstForward;
		
		ExtendedIterator<String[]> sit = iterateSequences(files);
		sit.progress(new ConsoleProgress(System.err), -1, r->"Processing reads").parallelized(nthreads, 1024, ei->ei.map(reads->{
			if (firstForward && reads[0]!=null) searchReads(reads[0],rnd,aho,false);
			if (firstReverse && reads[0]!=null) searchReads(SequenceUtils.getDnaReverseComplement(reads[0]),rnd,aho,false);
			if (firstForward && reads[1]!=null) searchReads(SequenceUtils.getDnaReverseComplement(reads[1]),rnd,aho,false);
			if (firstReverse && reads[1]!=null) searchReads(reads[1],rnd,aho,false);
			if (firstForward && reads[0]!=null) searchReads(reads[0],rnd,aho,true);
			if (firstReverse && reads[0]!=null) searchReads(SequenceUtils.getDnaReverseComplement(reads[0]),rnd,aho,true);
			if (firstForward && reads[1]!=null) searchReads(reads[1],rnd,aho,true);
			if (firstReverse && reads[1]!=null) searchReads(SequenceUtils.getDnaReverseComplement(reads[1]),rnd,aho,true);
			return null;
		})).drain();
		
		return null;
	}

	@SuppressWarnings("unchecked")
	private static ExtendedIterator<String[]> iterateSequences(String[] files) throws IOException {
		if (files.length==1) {
			boolean fastq = EI.lines(files[0]).first().startsWith("@");
			return EI.lines(files[0]).block(fastq?4:2,String.class).map(b->b[0].endsWith("/2")?new String[] {null,b[1]}:new String[] {b[1],null});
		}
		if (files.length!=2) throw new RuntimeException("Specify either one or two fasta or fastq files!");
		boolean fastq = EI.lines(files[0]).first().startsWith("@");
		return EI.lines(files[0]).block(fastq?4:2,String.class).
				fuse(String[].class, EI.lines(files[1]).block(fastq?4:2,String.class)).
				map(b->new String[] {b[0][1],b[1][1]});
	}


	private static boolean isDna(String seq) {
		return new Alphabet("ACGTNacgtn".toCharArray()).isValid(seq);
	}
	
	
	// frame 0-6: target, 6-12 decoy
	// frame odd: reverse
	private static String searchDna(String name, String seq, int frame, RandomNumbers rnd,
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, int minl, int maxl,
			MemoryIntervalTreeStorage<NameAnnotation> vars) {
	
		boolean transcr = ImmutableReferenceGenomicRegion.canParse(name);//name.contains(":");
		boolean decoy = frame>=6;
		boolean reverse = frame%2==1;
		ImmutableReferenceGenomicRegion<Object> parent = transcr?ImmutableReferenceGenomicRegion.parse(name):new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name, !reverse), new ArrayGenomicRegion(0,seq.length()));
		
		boolean extra = frame>=102 && frame<=103;
		boolean rnaseq = frame>=202 && frame<=207;
		
		boolean spliced = frame>=42 && frame<=43;
		if (spliced) {
			decoy=frame==43;
			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, 0, false);
			if (decoy && rnd!=null)
				aa = rnd.shuffle(aa.toString());
			else if (decoy)
				aa = StringUtils.reverse(aa);

			findSpliced(aa.toString(),aho,parent,decoy,minl,maxl);
			return null;
		}
		
		boolean substitution = frame>=92 && frame<=93;
		if (substitution) {
			decoy=frame==93;
			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, 0, false);
			if (decoy && rnd!=null)
				aa = rnd.shuffle(aa.toString());
			else if (decoy)
				aa = StringUtils.reverse(aa);

			findSubstitution(aa.toString(),aho,parent,decoy,maxl);
			return null;
		}
		
		boolean frameshift = frame>=84 && frame<=85;
		if (frameshift) {
			decoy=frame==85;
			int stop = SequenceUtils.translate(seq).indexOf('*');
			if (stop==-1) throw new RuntimeException("No Stop for \n>"+name+"\n"+seq);
			
			String[] triplets = new String[seq.length()/3];
			for (int i=0; i<triplets.length; i++)
				triplets[i] = seq.substring(i*3,i*3+3);
			
			if (decoy && rnd!=null){
				rnd.shuffle(triplets,0,stop);
				rnd.shuffle(triplets,stop+1,triplets.length);
			}
			else if (decoy){
				ArrayUtils.reverse(triplets, 0,stop);
				ArrayUtils.reverse(triplets, stop+1,triplets.length);
			}
			
			seq = EI.wrap(triplets).concat();

			findFrameshift(seq,stop*3,aho,parent,decoy,maxl);
			return null;
		}

		
		
		if (rnaseq) {
			frame = (frame-202)%3;
			decoy = (frame>=205);
			reverse = false;
		}
		else
			frame = (frame%6)/2;
		

		
		if (transcr && reverse) return null;

		// new: CharDag
		if (frame!=0) return null;
		
		if (reverse)
			seq = SequenceUtils.getDnaReverseComplement(seq);
		
		CharSequence sseq;
		if (decoy && rnd!=null)
			sseq = rnd.shuffle(seq.toString());
		else if (decoy)
			sseq = StringUtils.reverse(seq);
		else 
			sseq = seq;
		
		ImmutableReferenceGenomicRegion<Object> iparent = parent.toMutable().toStrandIndependent().toImmutable();
		int[] c = {0};
		CharDag dag = new CharDag(CharIterator.fromCharSequence(sseq),
					vars==null?EI.empty():vars.ei(iparent)
						.filter(v->iparent.contains(v))
						.unfold(r->{
							ArrayGenomicRegion reg = parent.induce(r.getRegion());
							String[] fromto = StringUtils.split(r.getData().getName(), '>');
							char[] from = parent.getReference().isMinus()?SequenceUtils.getDnaReverseComplement(fromto[0]).toCharArray():fromto[0].toCharArray();
							return EI.split(fromto[1], ',').map(ft->{
								char[] to = parent.getReference().isMinus()?SequenceUtils.getDnaReverseComplement(ft).toCharArray():ft.toCharArray();
								return new CharDag.SequenceVariant(reg.getStart(), from,to, r.getReference()+":"+r.getRegion().getStart()+" "+String.valueOf(from)+">"+ft).compact();
							});
						})
					);

		Trie<String> codeTrie = SequenceUtils.codeTrie;
		
		boolean udecoy = decoy;
		boolean ureverse = reverse;
		dag.traverse(i->i<3?new TranslateCharDagVisitor<>(aho.ahoCorasickVisitor(l->l*3),codeTrie,udecoy):null,(res, varit)->{
			int s;
			if (!udecoy)
				s = res.getStart();
			else
				s = (sseq.length()-res.getEnd());
			
			String varstring = EI.wrap(varit).concat(" ");
				
			if (ureverse) s = sseq.length()-s-res.getLength();
			if (transcr) {
				GenomicRegion reg = new ArrayGenomicRegion(s,s+Math.max(1, res.getLength()));
				if (parent.getRegion().getTotalLength()<=reg.getStop()) {
					System.out.println(res);
					throw new RuntimeException();
				}
				try {
				reg = parent.map(reg);
				} catch (IllegalArgumentException e) {
					e.printStackTrace();
					System.err.println(sseq);
					System.err.println(parent);
					System.err.println(reg);
					System.err.println(res);
					System.err.println(varstring);
				}
				if (reg.getNumParts()>1) 
					res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(udecoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()+","+varstring));
				// otherwise it's in the genome!
			} else if (extra) {
				res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(udecoy?"REV_"+name:name,!ureverse), new ArrayGenomicRegion(s,s+res.getLength()), res.getKey().toString()+",Extra"));
			} else if (rnaseq) {
				res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(udecoy?"REV_"+name:name,!ureverse), new ArrayGenomicRegion(s,s+res.getLength()), res.getKey().toString()+",RNA-seq"));
			} else {
				res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(udecoy?"REV_"+name:name,!ureverse), new ArrayGenomicRegion(s,s+Math.max(1, res.getLength())), res.getKey().toString()+","+varstring));
			}
			c[0]++;
		});
		
////			CharSequence aa = SequenceUtils.translate((!reverse?seq:SequenceUtils.getDnaReverseComplement(seq)).substring(frame));
//			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, frame, reverse);
//			if (decoy && rnd!=null)
//				aa = rnd.shuffle(aa.toString());
//			else if (decoy)
//				aa = StringUtils.reverse(aa);
//			
//			int c = 0;
//			for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(aa).loop()) {
//				int s;
//				if (!decoy)
//					s = res.getStart()*3+frame;
//				else
//					s = (aa.length()-res.getEnd())*3+frame;
//				
//				if (reverse) s = seq.length()-s-res.getLength()*3;
//				if (transcr) {
//					GenomicRegion reg = parent.map(new ArrayGenomicRegion(s,s+res.getLength()*3));
//					if (reg.getNumParts()>1) 
//						res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()+","));
//					// otherwise it's in the genome!
//				} else if (extra) {
//					res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,!reverse), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()+",Extra"));
//				} else if (rnaseq) {
//					res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,!reverse), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()+",RNA-seq"));
//				} else {
//					res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,!reverse), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()+","));
//				}
//			c++;
//			}
		if (transcr || rnaseq)
			return null;
		return "For sequence "+name+" reverse "+reverse+" decoy "+decoy+": "+c[0];
	}
	
	private static ImmutableReferenceGenomicRegion<String> readDummyRegion = new ImmutableReferenceGenomicRegion<String>(Chromosome.obtain("READ"),new ArrayGenomicRegion(0,1),",");
	private static ImmutableReferenceGenomicRegion<String> readDecoyDummyRegion = new ImmutableReferenceGenomicRegion<String>(Chromosome.obtain("REV_READ"),new ArrayGenomicRegion(0,1),",");
	private static void searchReads(String seq, RandomNumbers rnd,
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, boolean decoy) {
	
		
		
		CharSequence sseq;
		if (decoy && rnd!=null)
			sseq = rnd.shuffle(seq.toString());
		else if (decoy)
			sseq = StringUtils.reverse(seq);
		else 
			sseq = seq;
		
		int[] c = {0};
		CharDag dag = new CharDag(CharIterator.fromCharSequence(sseq),EI.empty());

		Trie<String> codeTrie = SequenceUtils.codeTrie;
		
		dag.traverse(i->i<3?new TranslateCharDagVisitor<>(aho.ahoCorasickVisitor(l->l*3),codeTrie,decoy):null,(res, varit)->{
			res.getValue().add(decoy?readDecoyDummyRegion:readDummyRegion); 
			c[0]++;
		});
	}
		
		
	private static String countDna(String name, String seq, int frame, 
			PeptideCounter counter, int minl, int maxl) {
	
		boolean transcr = name.contains(":");
		boolean decoy = frame>=6;
		boolean reverse = frame%2==1;
		ImmutableReferenceGenomicRegion<Object> parent = transcr?ImmutableReferenceGenomicRegion.parse(name):new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name, !reverse), new ArrayGenomicRegion(0,seq.length()));
		
		
		boolean extra = frame>=102 && frame<=103;
		boolean rnaseq = frame>=202 && frame<=207;
		
		boolean spliced = frame>=42 && frame<=43;
		if (spliced) {
			decoy=frame==43;
			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, 0, false);
			if (decoy) return null;

			countSpliced(aa.toString(),counter,parent,minl,maxl);
			return null;
		}
		
		boolean substitution = frame>=92 && frame<=93;
		if (substitution) {
			return null;
		}
		
		boolean frameshift = frame>=84 && frame<=85;
		if (frameshift) {
			return null;
		}

		
		
		if (rnaseq) {
			frame = (frame-202)%3;
			decoy = (frame>=205);
			reverse = false;
		}
		else
			frame = (frame%6)/2;
		
		if (decoy) return null;
		
		
		if (transcr && reverse) return null;

		
//			CharSequence aa = SequenceUtils.translate((!reverse?seq:SequenceUtils.getDnaReverseComplement(seq)).substring(frame));
			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, frame, reverse);
			int uframe = frame;
			boolean ureverse = reverse;
			int[] c = {0};
			counter.count(aa, (s,word)->{
				s = s*3+uframe;
				c[0]++;
				if (ureverse) s = seq.length()-s-word.length()*3;
				if (transcr) {
					GenomicRegion reg = parent.map(new ArrayGenomicRegion(s,s+word.length()*3));
					if (reg.getNumParts()>1) 
						return new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(parent.getReference().getName(),parent.getReference().getStrand()), reg, word.toString()+",");
					// otherwise it's in the genome!
				} else if (extra) {
					return new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name,!ureverse), new ArrayGenomicRegion(s,s+word.length()*3), word.toString()+",Extra");
				} else if (rnaseq) {
					return new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name,!ureverse), new ArrayGenomicRegion(s,s+word.length()*3), word.toString()+",RNA-seq");
				} else {
					return new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name,!ureverse), new ArrayGenomicRegion(s,s+word.length()*3), word.toString()+",");
				}
				return null;
			});
		if (transcr || rnaseq)
			return null;
		return "For sequence "+name+" reverse "+reverse+" decoy "+decoy+": "+c[0];
	}

	
	public static void findFrameshift(String seq, int stop, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho,
			ImmutableReferenceGenomicRegion<Object> parent, boolean decoy, int maxl) {
	
		for (int shift : new int[]{-2,-1,1,2}) {
			for (int s=3; s<stop; s+=3) {
				
				String before = seq.substring(Math.max(0, s-maxl*3),s);
				String after = seq.substring(s+shift,Math.min(seq.length(), s+shift+maxl*3));
				ArrayGenomicRegion bareg = new ArrayGenomicRegion(s-before.length(),s+shift+after.length());
				String aa = SequenceUtils.translate(before+after);
				if (aa.contains("*")) aa = aa.substring(0,aa.indexOf('*'));

				
				for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res : aho.iterateAhoCorasick(aa).loop()) {
					
					if (res.getStart()<before.length()/3 & res.getEnd()>before.length()/3)  {
						
						GenomicRegion reg = res.asRegion();
						reg = reg.pep2dna();
						reg = reg.extendBack(-1);
						reg = bareg.map(reg);
						reg = parent.map(reg);
						
						res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()+","+shift));
						
					}
					
				}
			}
		}
		
	}
	
	private static String[] aas = EI.wrap(new HashSet<String>(SequenceUtils.code.values())).filter(s->!s.equals("L")).toArray(new String[0]);
	
	public static void findSubstitution(String aa, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho,
			ImmutableReferenceGenomicRegion<Object> parent, boolean decoy, int maxl) {
	
		for (int p=0; p<aa.length(); p++) {
			
			String before = aa.substring(Math.max(0, p-maxl),p);
			String after = aa.substring(p+1,Math.min(aa.length(), p+1+maxl));
			ArrayGenomicRegion bareg = new ArrayGenomicRegion(p-before.length(),p+1+after.length());
			
			for (String subs : aas) {
				if (subs.equals(aa.substring(p,p+1)) || subs.equals("I")) continue;
				String saa = before+subs+after;
				
				for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res : aho.iterateAhoCorasick(saa).loop()) {
					
					if (res.getStart()<=before.length() & res.getEnd()>before.length())  {
						
						GenomicRegion reg = res.asRegion();
						reg = bareg.map(reg);
						reg = reg.pep2dna();
						reg = parent.map(reg);
						String subsstring = (before.length()-res.getStart())+":"+aa.substring(p,p+1)+">"+subs;
						res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()+","+subsstring));
						
					}
					
				}
			}
		}
		
	}
	
	private static void findSlipperyFrameshift(String seq, int stop, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho,
			ImmutableReferenceGenomicRegion<Object> parent, boolean decoy, int maxl) {
		
		Trie<Integer> slippery = new Trie<Integer>();
		// XXX YYY Z (X=A, G. U; Y=A, U; Z=A, C, U)
		for (String x : new String[]{"A","C","G","T"})
			for (String y : new String[]{"A","C","G","T"})
				for (String z : new String[]{"A","C","G","T"})
					slippery.put(x+x+x+y+y+y+z, -1);
		
		for (AhoCorasickResult<Integer> r : slippery.iterateAhoCorasick(seq).loop()) {
			if (r.getStart()%3==2 && r.getStart()<stop) {
				String before = seq.substring(Math.max(0, r.getEnd()-maxl*3),r.getEnd());
				String after = seq.substring(r.getStop(),Math.min(seq.length(), r.getStop()+maxl*3));
				ArrayGenomicRegion bareg = new ArrayGenomicRegion(r.getEnd()-before.length(),r.getStop()+after.length());
				String aa = SequenceUtils.translate(before+after);
				if (aa.contains("*")) aa = aa.substring(0,aa.indexOf('*'));

				for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res : aho.iterateAhoCorasick(aa).loop()) {
					
					if (res.getStart()<before.length()/3 & res.getEnd()>before.length()/3)  {
						
						GenomicRegion reg = res.asRegion();
						reg = reg.pep2dna();
						reg = reg.extendBack(-1);
						reg = bareg.map(reg);
						reg = parent.map(reg);
						
						res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()));
						
					}
					
				}
			}
		}
	}
	
	public static void findSpliced(String aa, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho,
			ImmutableReferenceGenomicRegion<Object> parent, boolean decoy, int minl, int maxl) {
		
		int inter = 25;
		SplitSequence split = new SplitSequence(aa, 0, 1, 2, maxl-1);
		
		for (int s=1; s<aa.length()+1; s++) {
			split.len1 = Math.min(s, maxl-1);
			split.start1 = s-split.len1;
//				for (int l=1; l<inter && s+l<aa.length()-1; l++) {
//					split.start2=s+l;
//					split.len2 = Math.min(maxl-1, aa.length()-s-l);
//					
//					findSpliced(aa, aho, parent, decoy, maxl);
//				}
			// in-order split
			for (int s2=s+1; s2-s<=inter && s2<aa.length()-1; s2++) {
				split.start2=s2;
				split.len2 = Math.min(maxl-1, aa.length()-s2);
				
				findSplit(split, aho, parent, decoy);
			}
			// reverse-order split is 1. difficult to handle, 2. difficult to store in the rgr and 3. useless anyways
			for (int s2=s-minl; s-s2<=maxl+inter && s2>=0; s2--) {
				split.start2=s2;
				split.len2 = Math.min(maxl-1, s-s2);
				
				findSplit(split, aho, parent, decoy);
			}
		}
		
	}
	
	
	public static void findSplit(SplitSequence split, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, ImmutableReferenceGenomicRegion<Object> parent, boolean decoy) {
		for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(split).loop()) {
			if (res.getStart()<split.len1 & res.getEnd()>split.len1)  {
				
				
				GenomicRegion pleft = new ArrayGenomicRegion(
						res.getStart()+split.start1,
						split.start1+split.len1
						);
			
				GenomicRegion pright = new ArrayGenomicRegion(
						split.start2,
						split.start2+res.getLength()-split.len1+res.getStart()
						);
			
				GenomicRegion preg = pleft.pep2dna().union(pright.pep2dna());
				
				if (decoy)
					preg = preg.reverse(parent.getRegion().getTotalLength());
				GenomicRegion reg = parent.map(preg);
				
				String extra = pleft.getTotalLength()+":"+split.getDistance(res.getStart(),res.getEnd())+":"+pright.getTotalLength()+":"+split.getLeftConsec(pleft,pright)+":"+split.getRightConsec(pleft,pright);
				res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()+","+extra));
				
			}
			
		}
	}

	
	public static void countSpliced(String aa, PeptideCounter counter,
			ImmutableReferenceGenomicRegion<Object> parent, int minl, int maxl) {
		
		int inter = 25;
		SplitSequence split = new SplitSequence(aa, 0, 1, 2, maxl-1);
		
		for (int s=1; s<aa.length()+1; s++) {
			split.len1 = Math.min(s, maxl-1);
			split.start1 = s-split.len1;
//				for (int l=1; l<inter && s+l<aa.length()-1; l++) {
//					split.start2=s+l;
//					split.len2 = Math.min(maxl-1, aa.length()-s-l);
//					
//					findSpliced(aa, aho, parent, decoy, maxl);
//				}
			// in-order split
			for (int s2=s+1; s2-s<=inter && s2<aa.length()-1; s2++) {
				split.start2=s2;
				split.len2 = Math.min(maxl-1, aa.length()-s2);
				
				countSplit(split, counter, parent);
			}
			// reverse-order split is 1. difficult to handle, 2. difficult to store in the rgr and 3. useless anyways
			for (int s2=s-minl; s-s2<=maxl+inter && s2>=0; s2--) {
				split.start2=s2;
				split.len2 = Math.min(maxl-1, s-s2);
				
				countSplit(split, counter, parent);
			}
		}
		
	}
	
	
	public static void countSplit(SplitSequence split, PeptideCounter counter, ImmutableReferenceGenomicRegion<Object> parent) {
		
		counter.count(split, (s,word)-> {
			if (s<split.len1 & s+word.length()>split.len1)  {
				
				
				GenomicRegion pleft = new ArrayGenomicRegion(
						s+split.start1,
						split.start1+split.len1
						);
			
				GenomicRegion pright = new ArrayGenomicRegion(
						split.start2,
						split.start2+word.length()-split.len1+s
						);
			
				GenomicRegion preg = pleft.pep2dna().union(pright.pep2dna());
				
				GenomicRegion reg = parent.map(preg);
				
				String extra = pleft.getTotalLength()+":"+split.getDistance(s,s+word.length())+":"+pright.getTotalLength()+":"+split.getLeftConsec(pleft,pright)+":"+split.getRightConsec(pleft,pright);
				return new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(parent.getReference().getName(),parent.getReference().getStrand()), reg, word.toString()+","+extra);
			}
			return null;
		});
		
	}


	private static class SplitSequence implements CharSequence {
		
		private String parent;
		private int start1;
		private int len1;
		private int start2;
		private int len2;

		public SplitSequence(String parent, int start1, int len1, int start2, int len2) {
			this.parent = parent;
			this.start1 = start1;
			this.len1 = len1;
			this.start2 = start2;
			this.len2 = len2;
		}

		public String getLeftConsec(GenomicRegion left, GenomicRegion right) {
			if (start2-left.getTotalLength()<0) return "";
			return parent.substring(start2-left.getTotalLength(),start2);
		}

		public String getRightConsec(GenomicRegion left, GenomicRegion right) {
			if (start1+len1+right.getTotalLength()>parent.length()) return "";
			return parent.substring(start1+len1,start1+len1+right.getTotalLength());
		}

		public int getDistance(int aastart, int aaend) {
			if (start2<start1)
				return -(aastart+ start1-(start2+aaend-len1));
			return start2-(start1+len1);
		}

		@Override
		public int length() {
			return len1+len2;
		}
		
		@Override
		public int hashCode() {
			return StringUtils.hashCode(this);
		}
		
		@Override
		public boolean equals(Object obj) {
			return StringUtils.equals(this, obj);
		}
		
		@Override
		public char charAt(int index) {
			return index<len1?parent.charAt(index+start1):parent.charAt(index-len1+start2);
		}
		
		@Override
		public String toString() {
			return StringUtils.toString(this);
		}

		@Override
		public CharSequence subSequence(int start, int end) {
			return toString().substring(start, end);
		}

	}

	private static String searchProtein(String name, String seq, int frame, RandomNumbers rnd,
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho) {
	
		boolean extra = frame>=102 && frame<=103;
		boolean decoy = frame%2==1;
		
		seq=seq.toUpperCase();
		
		
		if (decoy && rnd!=null)
			seq = rnd.shuffle(seq);
		else if (decoy)
			seq = StringUtils.reverse(seq).toString();
		
		String info = "";
		if (extra)
			info="Extra";
		else
			throw new RuntimeException("Cannot be!");
		
		for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(seq).loop()) {
			int s = res.getStart();
			if (decoy) s = seq.length()-s-res.getLength()*3;
			
			res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,true), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()+","+info));
		}
		return null;
	}
	private static String countProtein(String name, String seq, int frame,
			PeptideCounter counter) {
	
		boolean extra = frame>=102 && frame<=103;
		boolean decoy = frame%2==1;
		
		if (decoy) return null;
		
		
		String info = "Extra";
		if (!extra)
			throw new RuntimeException("Cannot be!");
		
		seq  = seq.replace('I', 'L');
		
		counter.count(seq, (s,word)->new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(name,true), new ArrayGenomicRegion(s,s+word.length()*3), word+","+info));
		return null;
	}

	
	private static boolean isDecoy(ImmutableReferenceGenomicRegion<String> hit) {
		return hit.getReference().getName().startsWith("REV");
	}
	
	
	public static class ModifiedAminoAcid {
		char aa;
		String delta;
		
		private void set(char aa, String delta) {
			this.aa = aa;
			this.delta = delta;
		}

		public static ModifiedAminoAcid parseSingle(String s) {

			Matcher m = Pattern.compile("^(.)\\(([+-.0-9]+)\\)$").matcher(s);
			if (!m.find()) throw new RuntimeException("Cannot parse "+s+" as modification!");
			ModifiedAminoAcid re = new ModifiedAminoAcid();
			re.aa = m.group(1).charAt(0);
			re.delta = m.group(2);
			return re;
		}
		
		public static int parse(String s, ModifiedAminoAcid[] buff) {

			int len = 0;
			int p = 0;
			for (int sepIndex=s.indexOf('('); sepIndex>=0; sepIndex = s.indexOf('(',p)) {
				for (int i=p; i<sepIndex-1; i++)
					buff[len++].set(s.charAt(i),null);
				buff[len++].set(s.charAt(sepIndex-1),s.substring(sepIndex+1, s.indexOf(')',p)));
				
				p = s.indexOf(')',p)+1;
			}
			for (int i=p; i<s.length(); i++)
				buff[len++].set(s.charAt(i),null);
			buff[len].set((char)0,null);
			return len;
		}

		public char getAlternative(ModifiedAminoAcid[] buff, int ind) {
			if (aa=='L') return 'I';
			if (aa=='I') return 'L';
			if (ind+1<buff.length && buff[ind+1].aa=='G' && aa=='D')
				return 'N';
			if (ind==0 && aa=='Q' && "-17.03".equals(delta))
				return 'E';
			return 0;
		}
		@Override
		public String toString() {
			return aa+"("+delta+")";
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + aa;
			result = prime * result + ((delta == null) ? 0 : delta.hashCode());
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
			ModifiedAminoAcid other = (ModifiedAminoAcid) obj;
			if (aa != other.aa)
				return false;
			if (delta == null) {
				if (other.delta != null)
					return false;
			} else if (!delta.equals(other.delta))
				return false;
			return true;
		}
		
	}
	
	private static boolean areEqualUptoIsobar(String a, String b) {
		if (a.length()!=b.length()) return false;
		a = a.replace('I', 'L');
		b = b.replace('I', 'L');
		if (a.startsWith("Q")) a = "E"+a.substring(1);
		if (b.startsWith("Q")) b = "E"+b.substring(1);
		a = a.replace("NG","DG");
		b = b.replace("NG","DG");
		return a.equals(b);
	}




	public static class FindGenomicPeptidesParameterSet extends GediParameterSet {

		public GediParameter<Integer> deltaFirst = new GediParameter<Integer>(this,"df", "The maximal score difference to the best sequence (must not have an identified location).", false, new IntParameterType(), 15);
		public GediParameter<Integer> deltaNext = new GediParameter<Integer>(this,"dn", "The minimal score difference to the next best sequence (with identified location)", false, new IntParameterType(), 16);
		public GediParameter<Integer> minlen = new GediParameter<Integer>(this,"minlen", "The minimal length of a peptide to consider", false, new IntParameterType(), 8);
		public GediParameter<Integer> maxlen = new GediParameter<Integer>(this,"maxlen", "The maximal length of a peptide to consider", false, new IntParameterType(), 22);
		public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
		public GediParameter<String> input = new GediParameter<String>(this,"in", "Peaks output csv", true, new StringParameterType());
		public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "Additional prefix", true, new StringParameterType(),"");
		public GediParameter<Genomic> genomic = new GediParameter<Genomic>(this,"g", "Genomic name", true, new GenomicParameterType());
		public GediParameter<String> tis = new GediParameter<String>(this,"tis", "Start codons to consider (prioritized)", true, new StringParameterType(),"AUG,CUG,ACG,GUG,AUC");

		public GediParameter<String> extra = new GediParameter<String>(this,"extra", "Fasta file containing extra sequences (DNA will be 6-frame translated!)!", false, new StringParameterType(),true);
		public GediParameter<String> rnaseq = new GediParameter<String>(this,"rnaseq", "Fasta file containing extra transcript sequences (will be 3-frame translated!)!", false, new StringParameterType(),true);
		public GediParameter<String> reads = new GediParameter<String>(this,"reads", "Fastq file(s) or Fasta file with RNA-seq reads (will be 3-frame translated!)!", true, new StringParameterType(),true);
		public GediParameter<Strandness> readsStrandness = new GediParameter<Strandness>(this,"strandness", "Strandness of the reads (cannot be automatic!)", false, new EnumParameterType<>(Strandness.class),true);
		
		public GediParameter<File> readSeqs = new GediParameter<File>(this,"${in}${prefix}.readseq", "Folder of identified reads sequences", false, new FileParameterType());
		public GediParameter<File> readSeqTab = new GediParameter<File>(this,"${in}${prefix}.readseq.tsv.gz", "Table of identified reads sequences", false, new FileParameterType());

		
		public GediParameter<String> hla = new GediParameter<String>(this,"hla", "File containing HLA allels (each line); could also be in ${input%.csv}.hla, in which case you do not have to specify this parameter!", true, new StringParameterType(),true);
		public GediParameter<String> variants = new GediParameter<String>(this,"var", "File containing variants (each line); could also be in ${input%.csv}.var, in which case you do not have to specify this parameter!", true, new StringParameterType(),true);
		public GediParameter<String> netmhc = new GediParameter<String>(this,"netmhc", "Command to call netmhc", true, new StringParameterType(),"netMHCpan");
		
		public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "Only search against chromosome 22", false, new BooleanParameterType());
		public GediParameter<Boolean> rnd = new GediParameter<Boolean>(this,"rnd", "Randomize decoys", false, new BooleanParameterType());

		public GediParameter<Boolean> noloc = new GediParameter<Boolean>(this,"noloc", "Output also hits without location into .annotated.csv.gz", false, new BooleanParameterType());
		public GediParameter<Boolean> all = new GediParameter<Boolean>(this,"all", "All possible annotations", false, new BooleanParameterType());
		public GediParameter<Boolean> hydro = new GediParameter<Boolean>(this,"hydro", "Query SRRCalc to compute hydrophobicities", false, new BooleanParameterType());
		public GediParameter<String> cat = new GediParameter<String>(this,"cat", "Which categories to search for (either a file, or directly comma-separated, order is important!), all="+StringUtils.concat(',',Category.values())+" default="+StringUtils.concat(',', defaultCategories), false, new StringParameterType(),"default");
		public GediParameter<String> dmod = new GediParameter<String>(this,"dmod", "Specify modifications that make sequences distinct (in the same format as in the Peptide column, separate multiple modifications by comma, e.g. -dmod S(+79.97),T(+79.97),Y(+79.97).", false, new StringParameterType(),"");
				
		public GediParameter<File> countOut = new GediParameter<File>(this,"${anchor}.count.tsv", "Output table for peptide count statistics", false, new FileParameterType());
		public GediParameter<String> anchors = new GediParameter<String>(this,"anchor", "Anchor residues for peptide count statistics", false, new StringParameterType());

		
		public GediParameter<File> annotatedPeaksOut = new GediParameter<File>(this,"${in}${prefix}.annotated.csv.gz", "Peaks file annotated with matched sequence, location and type (filtered by df and pep length)", false, new FileParameterType());
		public GediParameter<File> nextFilteredOut = new GediParameter<File>(this,"${in}${prefix}.deltaNext.csv.gz", "Peaks file annotated by delta next filter; filtered for best hit per Sequence", false, new FileParameterType());
		public GediParameter<Boolean> writeUnidentified = new GediParameter<Boolean>(this,"unidentified", "Write a file with all top hits (also the unidentified ones!)", false, new BooleanParameterType());
		public GediParameter<Boolean> writeNextFiltered = new GediParameter<Boolean>(this,"deltaNextStats", "Write a file with reports for delta next", false, new BooleanParameterType());
		
		public GediParameter<File> unidentifiedPeaksFOut = new GediParameter<File>(this,"${in}${prefix}.pep.unidentified.csv", "Peptide list of unannotated spectra", false, new FileParameterType());
		
		public GediParameter<File> pepGenomeOut = new GediParameter<File>(this,"${in}${prefix}.pep.tmp", "Peaks file filtered for best hit per Sequence (filtered by dn)", false, new FileParameterType()).setRemoveFile(true);

		public GediParameter<File> fdrOut = new GediParameter<File>(this,"${in}${prefix}.fdr.csv", "FDR statistics", false, new FileParameterType());
		public GediParameter<File> fdrPlot = new GediParameter<File>(this,"${in}${prefix}.fdr.pdf", "FDR statistics plot", false, new FileParameterType());

		public GediParameter<File> pepGenomeFOut = new GediParameter<File>(this,"${in}${prefix}.pep.csv.gz", "Peaks file filtered for best hit per Sequence (filtered by dn)", false, new FileParameterType());
		public GediParameter<File> mhcOut = new GediParameter<File>(this,"${in}${prefix}.pep.mhc.csv.gz", "netMHC predictions", false, new FileParameterType());
		public GediParameter<File> rtOut = new GediParameter<File>(this,"${in}${prefix}.pep.rt.csv.gz", "SRRCalc", false, new FileParameterType());
		
		public GediParameter<File> pepGenomeAOut = new GediParameter<File>(this,"${in}${prefix}.pep.annotated.csv.gz", "Annotated peaks file (Orfs, NetMHC binding)", false, new FileParameterType());
		
		public GediParameter<Double> pepBedOutThreshold = new GediParameter<Double>(this,"qbed", "Q value threshold for bed output", false, new DoubleParameterType(), 0.1);
		public GediParameter<Double> pepBedOutMHCThreshold = new GediParameter<Double>(this,"mhcbed", "netMHC precition value threshold for bed output", false, new DoubleParameterType(), 0.5);
		public GediParameter<File> pepBedOut = new GediParameter<File>(this,"${in}${prefix}.pep.bed", "Bed file of all identified peptides", false, new FileParameterType());

		public GediParameter<File> paramFile = new GediParameter<File>(this,"${in}${prefix}.param", "File containing the parameters used to call Prism", false, new FileParameterType());
		public GediParameter<File> runtimeFile = new GediParameter<File>(this,"${in}${prefix}.runtime", "File containing the runtime information", false, new FileParameterType());
		
	}

	
}
