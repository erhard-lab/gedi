package executables;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.logging.Level;

import cern.colt.bitvector.BitVector;
import gedi.core.data.annotation.Transcript;
import gedi.core.genomic.Genomic;
import gedi.core.reference.Chromosome;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.RunUtils;
import gedi.util.SequenceUtils;
import gedi.util.StringUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.tree.Trie;
import gedi.util.datastructure.tree.Trie.AhoCorasickResult;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineIterator;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.fasta.DefaultFastaHeaderParser;
import gedi.util.io.text.fasta.FastaEntry;
import gedi.util.io.text.fasta.FastaFile;
import gedi.util.io.text.fasta.FastaHeaderParser;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.mutable.MutableInteger;
import gedi.util.mutable.MutableMonad;
import gedi.util.mutable.MutablePair;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;
import gedi.util.r.RRunner;
import gedi.util.sequence.Alphabet;


public class FindGenomicPeptides2 {
	
	

	public static void main(String[] args) throws IOException {
		
		FindGenomicPeptidesParameterSet params = new FindGenomicPeptidesParameterSet();
		GediProgram pipeline = GediProgram.create("FindGenomicPeptides",
				new FindGenomicPeptidesProgram(params),
				new ToPeptideListProgram(params, params.pepCdsOut, LocationAnnotation.CDS,LocationAnnotation.PeptideSpliced),
				new ToPeptideListProgram(params, params.pepTranscriptomeOut, LocationAnnotation.CDS,LocationAnnotation.OffFrame, LocationAnnotation.ncRNA, LocationAnnotation.UTR3,LocationAnnotation.UTR5,LocationAnnotation.PeptideSpliced),
				new ToPeptideListProgram(params, params.pepGenomeOut, LocationAnnotation.values()),
				new PrepareFdrProgram(params),
				new ComputeFdrProgram(params, params.pepCdsOut, params.pepCdsFOut, "CDS", LocationAnnotation.CDS,LocationAnnotation.PeptideSpliced),
				new ComputeFdrProgram(params, params.pepTranscriptomeOut, params.pepTranscriptomeFOut, "Transcriptome", LocationAnnotation.CDS,LocationAnnotation.OffFrame, LocationAnnotation.ncRNA, LocationAnnotation.UTR3,LocationAnnotation.UTR5,LocationAnnotation.PeptideSpliced),
				new ComputeFdrProgram(params, params.pepGenomeOut, params.pepGenomeFOut, "Genome", LocationAnnotation.values()),
				new AnnotateOrfAndMhcBindingProgram(params, params.pepCdsFOut, params.pepCdsAOut, "CDS"),
				new AnnotateOrfAndMhcBindingProgram(params, params.pepTranscriptomeFOut, params.pepTranscriptomeAOut, "Transcriptome"),
				new AnnotateOrfAndMhcBindingProgram(params, params.pepGenomeFOut, params.pepGenomeAOut, "Genome"),
				new PlotProgram(params)
//				new ProcessUnidentifiedProgram(params)
				);
		GediProgram.run(pipeline, null, null, new CommandLineHandler("FindGenomicPeptides","FindGenomicPeptides identifies a list of peptides in the genomic 6-frame translation.",args));

	}
	
	public static class PeptideOrf implements Comparable<PeptideOrf>{
		String dna;
		int[] startPrioPos;
		GenomicRegion pepLocInDna;
		public PeptideOrf(String dna, GenomicRegion pepLocInDna, int[] startPrioPos) {
			this.dna = dna;
			this.startPrioPos = startPrioPos;
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

		public String getTis() {
			return startPrioPos==null?"-":dna.substring(startPrioPos[1], startPrioPos[1]+3);
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
	
	
	public static List<ImmutableReferenceGenomicRegion<PeptideOrf>> getOrfs(ImmutableReferenceGenomicRegion<?> pep, Genomic g, Trie<Integer> starts) {
		
		if (pep.getReference().getName().startsWith("REV_"))
			pep = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(pep.getReference().getName().substring(4),pep.getReference().getStrand()), pep.getRegion());
		
		ImmutableReferenceGenomicRegion<?> upep = pep;
		
		ArrayList<ImmutableReferenceGenomicRegion<PeptideOrf>> re = new ArrayList<>();
		int ex = 99;
		int lu = -1;
		int ld = -1;
		
		do  {
			int uex = ex;
			String seq = SequenceUtils.translate(g.getSequence(pep.toMutable().transformRegion(r->r.extendFront(uex).extendBack(uex).intersect(new ArrayGenomicRegion(0,g.getLength(upep.getReference().getName()))))));
			int up = seq.substring(0, ex/3).lastIndexOf('*');
			ld = seq.indexOf('*', ex/3)*3+3-ex-pep.getRegion().getTotalLength();
			lu=up<0?-1:(ex/3-up)*3-3;
			ex=ex*3;
		} while ((lu<0 || ld<0) && pep.getRegion().getStart()-ex/3>=0 && pep.getRegion().getEnd()+ex/3<g.getLength(pep.getReference().getName()));
		
		if (lu>=0 && ld>=0) {
			MutableReferenceGenomicRegion<?> orf = pep.toMutable().extendRegion(lu,ld);
			String seq = g.getSequence(orf).toString();
			
			re.add(new ImmutableReferenceGenomicRegion<>(orf.getReference(), orf.getRegion(), new PeptideOrf(seq,new ArrayGenomicRegion(lu,orf.getRegion().getTotalLength()-ld),findStart(seq,lu,starts))));
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
				
				re.add(new ImmutableReferenceGenomicRegion<>(t.getReference(), t.toMutable().extendRegion(-offset, 0).map(xorf), new PeptideOrf(seq,new ArrayGenomicRegion(lu,xorf.getTotalLength()-ld),findStart(seq,lu,starts))));
			}
		}
		
		return re;
	}
	
	/**
	 * [start codon prio,start pos]
	 * @param seq
	 * @param pep
	 * @param starts
	 * @return
	 */
	private static int[] findStart(String seq, int pep, Trie<Integer> starts) {
		ArrayList<int[]> list = starts.iterateAhoCorasick(seq)
				.filter(res->res.getStart()%3==0 & res.getStart()<=pep)
				.map(res->new int[]{res.getValue(),res.getStart()})
				.list();
		if (list.size()==0) return null;
		list.sort((a,b)->{
			int re = Integer.compare(a[0], b[0]); // first by start codon priority
			if (re==0)
				re = Integer.compare(b[1], a[1]); // then by shorter ORF
			return re;
		});
		return list.get(0);
	}
	
	public static class Affinity implements Comparable<Affinity> {
		String allele;
		String sequence;
		double perrank;
		
		public Affinity(String allele, String sequence, double perrank) {
			this.allele = allele;
			this.sequence = sequence;
			this.perrank = perrank;
		}

		@Override
		public int compareTo(Affinity o) {
			return Double.compare(perrank, o.perrank);
		}
	}

	public static enum LocationAnnotation {
		CDS,UTR5,UTR3,OffFrame,ncRNA,Intronic,Intergenic,PeptideSpliced;
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
	
	public static ExtendedIterator<String[][]> iteratePeaksBlocks(String file, MutableMonad<HeaderLine> header, Consumer<String> first, Predicate<String[]> filterHit, boolean removeDoubleSeqAndResort) throws IOException {
		LineIterator lit = EI.lines(file);
		String h = lit.next();
		if (header!=null)
			header.Item = new HeaderLine(h,',');
		if (first!=null)
			first.accept(h);
		
		int alc = header.Item.get("ALC (%)");
		int fraction = header.Item.get("Fraction");
		int scan = header.Item.get("Scan");
		int peptide = header.Item.get("Peptide");
		
		
		Comparator<String[]> compFields = (a,b)->{
			int re = Integer.compare(Integer.parseInt(a[fraction]),Integer.parseInt(b[fraction]));
			if (re>0) return -1;
			if (re==0)
				re = Integer.compare(Integer.parseInt(a[scan]),Integer.parseInt(b[scan]));
			return re;
		};
		BitVector bv = new BitVector(64);
		HashSet<String> peps = new HashSet<>();
		
		
		ExtendedIterator<String[][]> re = lit.map(a->StringUtils.split(a, ','))
				.iff(filterHit!=null, ei->ei.filter(filterHit))
				.multiplex(compFields,String[].class);
		
		if (removeDoubleSeqAndResort)
			re = re.map(f->{
						Arrays.sort(f,(a,b)->Integer.compare(Integer.parseInt(b[alc]), Integer.parseInt(a[alc])));
						// remove double sequences (e.g. due to hcd and hcd/etd)
						bv.clear();
						if (bv.size()<f.length)
							bv.setSize(f.length);
						peps.clear();
						for (int i=0; i<f.length; i++) 
							bv.putQuick(i,peps.add(f[i][peptide]));
						if (bv.cardinality()<f.length)
							f = ArrayUtils.restrict(f, bv);
						return f;
					});
		return re;
	}
	
	public static ExtendedIterator<Affinity> iterateNetMHC(String netmhc, String hla, int l, String path) {
		return new LineIterator(RunUtils.output(netmhc,"-a",hla,"-inptype","1","-l",l+"","-f",path))
			.filter(s->s.length()>0 && !s.startsWith("#") && !s.startsWith("-") && !s.startsWith("Number"))
			.map(s->StringUtils.trim(s))
			.filter(s->!s.toLowerCase().startsWith("pos") && !s.toLowerCase().startsWith("protein") && !s.toLowerCase().contains("distance"))
			.map(s->s.split("\\s+"))
			.map(a->new Affinity(a[1],a[2],Double.parseDouble(a[12])));
	}

	
	
	public static class ProcessUnidentifiedProgram extends GediProgram {


		public ProcessUnidentifiedProgram(FindGenomicPeptidesParameterSet params) {
			
			addInput(params.unidentifiedPeaksOut);
			addInput(params.plot);
			addInput(params.hla);
			addInput(params.netmhc);
			addInput(params.fdrOut);
			
			addInput(params.input);
			
			// just such that the fdrs are computed
			addInput(params.pepCdsFOut);
			addInput(params.pepTranscriptomeFOut);
			addInput(params.pepGenomeFOut);
			
			addOutput(params.unidentifiedPeaksFOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			boolean plot = getBooleanParameter(1);
			String hla = getParameter(2);
			String netmhc = getParameter(3);
			File fdr = getParameter(4);
			
			String prefix = getParameter(5);
			
			String input = csv.getPath();
			
			
			HashMap<Integer,double[][]> fdrs = new HashMap<>();// len->[score]->[q1,q2]
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			for (String[] a : EI.lines(fdr.getPath())
				.skip(1, s->header.Item = new HeaderLine(s,','))
				.map(s->StringUtils.split(s, ','))
				.filter(s->s[header.Item.get("Mode")].equals("Genome")).loop()) {
				
				int len = Integer.parseInt(a[header.Item.get("Peptide length")]);
				int alc = Integer.parseInt(a[header.Item.get("ALC")]);
				double q1 = Double.parseDouble(a[header.Item.get("Q1")]);
				double q2 = Double.parseDouble(a[header.Item.get("Q2")]);
				
				fdrs.computeIfAbsent(len, x->new double[101][])[alc] = new double[]{q1,q2};
			}
			
			
			if (hla==null && new File(prefix+".hla").exists()) 
				hla = prefix+".hla";
			
			if (hla==null && new File(StringUtils.removeFooter(prefix,".csv")+".hla").exists()) 
				hla = StringUtils.removeFooter(prefix,".csv")+".hla";
			
			if (hla!=null) 
				hla = EI.lines(hla, "#").filter(s->s.length()>0).concat(",");
			
			
			HashMap<String,Affinity> aff = null;
			
			if (hla!=null) {
				try {
					if (new ProcessBuilder(netmhc,"-h").start().waitFor()==0) {
						context.getLog().info("netMHC works, will annotate binding affinities! (unidentified)");
						aff=new HashMap<>();
					}
					else
						context.getLog().info("netMHC does not work, will not annotate binding affinities! (unidentified)");
				} catch (IOException e) {
					context.getLog().log(Level.WARNING,"netMHC does not work, will not annotate binding affinities! (unidentified)",e);
				}
			} else
				context.getLog().info("No HLA alleles given, will not annotate binding affinities! (unidentified)");
			
			
			
			
			HashMap<String,int[]> best = new HashMap<>();
			
			HashMap<Integer,LineOrientedFile> pepout = new HashMap<>();
			for (String[] a : iteratePeaksLines(input, header,null).loop()) {
				String pep = removeMod(a[header.Item.get("Peptide")]);
				
				int alc = Integer.parseInt(a[header.Item.get("ALC (%)")]);
				int[] pp = best.computeIfAbsent(pep, x->new int[2]);
				pp[0] = Math.max(pp[0],alc);
				pp[1]++;
				
				if (aff!=null && pep.length()<=15 && pep.length()>=8) // outside of this range, netMHC does not work!
					pepout.computeIfAbsent(pep.length(),x->{
						try {
							LineOrientedFile re = new LineOrientedFile(File.createTempFile("FindGenomicPeptides.unidentified."+pep.length(), ".peplist").getPath());
							re.startWriting();
							return re;
						} catch (IOException e) {
							throw new RuntimeException("Could not write to temporary file!",e);
						}
					}).writeLine(pep+"\n"+StringUtils.reverse(pep)); // also compute decoys!
			}
			for (LineOrientedFile lo : pepout.values()) lo.finishWriting();
		
			if (aff!=null)
				for (Integer l : pepout.keySet()) {
					context.getLog().info("netMHC predictions for length "+l+" (unidentified)");
					
					for (Affinity a : iterateNetMHC(netmhc,hla,l,pepout.get(l).getPath()).loop()) {
						try{
							aff.merge(a.sequence, a, (x,y)->x.perrank<y.perrank?x:y);
						} catch (NumberFormatException e) {
							context.getLog().warning("Unexpected line: "+StringUtils.toString(a));
						}
					}
				}
			for (LineOrientedFile lo : pepout.values()) lo.delete();
			
			
			
			context.getLog().info("Annotate list for unidentified");

			HashMap<String, Affinity> uaff = aff;
			LineWriter out = getOutputWriter(0);
			for (String[] a : iteratePeaksLines(input, header,
													h->out.writeLine2(h+",Q1,Q2,spectra"+(uaff!=null?",HLA allele,netMHC % rank,Decoy HLA allele,Decoy netMHC % rank":"")) 
													).loop()) {
				
				String pep = removeMod(a[header.Item.get("Peptide")]);
				
				String affout = "";
				if (aff!=null) {
					Affinity af = aff.get(pep);
					if (af==null)
						affout=",-,100";
					else
						affout=String.format(",%s,%.2f", af.allele,af.perrank);
					
					af = aff.get(StringUtils.reverse(pep).toString());
					if (af==null)
						affout+=",-,100";
					else
						affout+=String.format(",%s,%.2f", af.allele,af.perrank);
				}
				
				
				boolean thereIsaBetterSpectrum = false;
				int[] bestHit = best.get(pep);
				if (bestHit==null || bestHit[0]!=Integer.parseInt(a[header.Item.get("ALC (%)")]))
					thereIsaBetterSpectrum = true;
				
				// remove it from the map (important if there is another spectrum with equal score)!
				if (!thereIsaBetterSpectrum) 
					best.remove(pep);
				
				
				if (!thereIsaBetterSpectrum) {
					
					int alc = Integer.parseInt(a[header.Item.get("ALC (%)")]);
					if (fdrs.containsKey(pep.length())) {
						double[] q = fdrs.get(pep.length())[alc];
						
						out.writef("%s,%.3g,%.3g,%d%s\n",StringUtils.concat(",", a),q[0],q[1],bestHit[1],affout);
					}
				}
				
				
			}
			
			
			out.close();

			if (plot) {
				try {
					context.getLog().info("Running R script for plotting");
					RRunner r = new RRunner(prefix+".mhc.unidentified.R");
					r.set("file",getOutputFile(0).getPath());
					r.addSource(getClass().getResourceAsStream("/resources/netMHC.R"));
					r.run(true);
				} catch (Throwable e) {
					context.getLog().log(Level.SEVERE, "Could not plot!", e);
				}
			}
			return null;
		}
		
	}
	
	
	public static ReentrantLock fdrLock = new ReentrantLock();
	
	public static class PlotProgram extends GediProgram {

		public PlotProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.input);
			addInput(params.pepCdsFOut);
			addInput(params.pepTranscriptomeFOut);
			addInput(params.pepGenomeFOut);
			
			setRunFlag(params.plot);
			
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
	
	
	public static class AnnotateOrfAndMhcBindingProgram extends GediProgram {

		private String mode;

		public AnnotateOrfAndMhcBindingProgram(FindGenomicPeptidesParameterSet params, GediParameter<File> in, GediParameter<File> out, String mode) {
			this.mode = mode;
			
			addInput(in);
			addInput(params.plot);
			addInput(params.tis);
			addInput(params.seq);
			addInput(params.hla);
			addInput(params.netmhc);
			addInput(params.spliced);
			
			addInput(params.input);
			
			addOutput(out);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			boolean plot = getBooleanParameter(1);
			String tis = getParameter(2);
			String[] fasta = getParameters(3).toArray(new String[0]);
			String hla = getParameter(4);
			String netmhc = getParameter(5);
			boolean spliced = getParameter(6);
			
			String prefix = getParameter(7);
			
			String input = csv.getPath();
			
			Trie<Integer> starts = new Trie<>();
			for (String t : EI.split(tis, ',').map(a->a.toUpperCase().replace('U', 'T')).loop())
				starts.put(t, starts.size());
			
			Genomic genomic;
			
			if (Genomic.check(fasta[0])!=null) {
				genomic = Genomic.get(fasta);
				fasta = Genomic.get(fasta).getGenomicFastaFiles().toArray(String.class);
				context.getLog().info("Recognized genomic input. Will annotate ORFs! ("+mode+")");
			} else {
				genomic = null;
				context.getLog().info("Recognized fasta. Won't annotate ORFs! ("+mode+")");
			}
			
			if (hla==null && new File(prefix+".hla").exists()) 
				hla = prefix+".hla";
			
			if (hla==null && new File(StringUtils.removeFooter(prefix,".csv")+".hla").exists()) 
				hla = StringUtils.removeFooter(prefix,".csv")+".hla";
			
			if (hla!=null) 
				hla = EI.lines(hla, "#").filter(s->s.length()>0).concat(",");
			
			
			HashMap<String,Affinity> aff = null;
			
			if (hla!=null) {
				try {
					if (new ProcessBuilder(netmhc,"-h").start().waitFor()==0) {
						context.getLog().info("netMHC works, will annotate binding affinities! ("+mode+")");
					
						// netMHC predictions
						HashMap<Integer,LineOrientedFile> pepout = new HashMap<>();
						MutableMonad<HeaderLine> header = new MutableMonad<>(); 
						for (String[] a : iteratePeaksLines(input, header,null).loop()) {
							String pep = removeMod(a[header.Item.get("Peptide")]);
							if (pep.length()<=15 && pep.length()>=8) // outside of this range, netMHC does not work!
								pepout.computeIfAbsent(pep.length(),x->{
									try {
										LineOrientedFile re = new LineOrientedFile(File.createTempFile("FindGenomicPeptides."+mode+"."+pep.length(), ".peplist").getPath());
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
							context.getLog().info("netMHC predictions for length "+l+" ("+mode+")");
							for (Affinity a : iterateNetMHC(netmhc,hla,l,pepout.get(l).getPath()).loop()) {
								try{
									aff.merge(a.sequence, a, (x,y)->x.perrank<y.perrank?x:y);
								} catch (NumberFormatException e) {
									context.getLog().warning("Unexpected line: "+StringUtils.toString(a));
								}
							}

						}
						
					}
					else
						context.getLog().info("netMHC does not work, will not annotate binding affinities! ("+mode+")");
				} catch (IOException e) {
					context.getLog().log(Level.WARNING,"netMHC does not work, will not annotate binding affinities! ("+mode+")",e);
				}
			} else
				context.getLog().info("No HLA alleles given, will not annotate binding affinities! ("+mode+")");
			
			
			
			
			
			RandomNumbers rnd = new RandomNumbers(42);
			
			context.getLog().info("Annotate list for "+mode);
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 

			HashMap<String, Affinity> uaff = aff;
			LineWriter out = getOutputWriter(0);
			for (String[] a : iteratePeaksLines(input, header,
													h->out.writeLine2(h+(genomic!=null?",ORF TIS,ORF length,ORF nterm,ORF cterm":"")+(uaff!=null?",HLA allele,netMHC % rank":"")+(spliced?",Left,Between,Right":"")) 
													).loop()) {
				
				String orf = "";
				String affout = "";
				String splic = "";
				if (genomic!=null) {
					String[] locs = StringUtils.split(a[header.Item.get("Location")],';');
					
					IntArrayList best = new IntArrayList();
					List<ImmutableReferenceGenomicRegion<PeptideOrf>> bestOrfs = new ArrayList<>();
					
					for (int i=0; i<locs.length; i++) {
						List<ImmutableReferenceGenomicRegion<PeptideOrf>> orfs = getOrfs(ImmutableReferenceGenomicRegion.parse(locs[i]), genomic, starts);
						orfs.sort((x,y)->x.getData().compareTo(y.getData()));
						
						if (orfs.size()>0) {
							int cmp = -1;
							if (bestOrfs.isEmpty() || (cmp=orfs.get(0).getData().compareTo(bestOrfs.get(0).getData()))<=0) {
								if (cmp<0) {
									best.clear();
									bestOrfs.clear();
								}
								best.add(i);
								bestOrfs.add(orfs.get(0));  // add just one of the potentially multiple best orfs for this location
							}
						}
					}

					String[] ann = StringUtils.split(a[header.Item.get("Annotation")],';');
					String[] seqs = StringUtils.split(a[header.Item.get("Sequence")],';');
					
					int select = 0;
					if (best.size()>1) 
						select=rnd.getUnif(0, best.size());
					if (best.size()==0) 
						select=rnd.getUnif(0, ann.length);
					
					
					a[header.Item.get("Sequence")]=seqs[select];
					a[header.Item.get("Location")]=locs[select];
					a[header.Item.get("Annotation")]=ann[select];
					a[header.Item.get("UniqueAnnotation")]=ann[select];
					
					if (best.size()==0) {
						orf=",-,-1,-1,-1";
					} else {
						PeptideOrf porf = bestOrfs.get(select).getData();
						orf=String.format(",%s,%d,%d,%d", porf.getTis(),porf.getLength(),porf.getNterm(),porf.getCterm());
					}
				}
				
				if (aff!=null) {
					String pep = removeMod(a[header.Item.get("Peptide")]);
					Affinity af = aff.get(pep);
					if (af==null)
						affout=",-,100";
					else
						affout=String.format(",%s,%.2f", af.allele,af.perrank);
					
					
				}
				
				if (spliced) {
					int left=0,right=0,between=0;
					
					if (a[header.Item.get("Category")].equals(LocationAnnotation.PeptideSpliced.name())) {
						ImmutableReferenceGenomicRegion<Object> loc = ImmutableReferenceGenomicRegion.parse(a[header.Item.get("Location")].replace("REV_", ""));
						
						ImmutableReferenceGenomicRegion<Transcript> trans = genomic.getTranscripts().ei(loc).filter(t->t.getRegion().contains(loc.getRegion()) && t.induce(loc.getRegion()).getNumParts()==2).first();
						
						if (trans==null) throw new RuntimeException("No transcript found for "+loc);
						
						ArrayGenomicRegion indu = trans.induce(loc.getRegion());
						left = indu.getPart(0).length()/3;
						right = indu.getPart(1).length()/3;
						between = indu.getIntronLength(0)/3;
					}
					
					splic=String.format(",%d,%d,%d",left,between,right);
				}
				
				out.writef("%s%s%s%s\n",StringUtils.concat(",", a),orf,affout,splic);
			}
			
			
			out.close();

			if (plot) {
				try {
					context.getLog().info("Running R script for plotting");
					RRunner r = new RRunner(prefix+".mhc."+mode+".R");
					r.set("file",getOutputFile(0).getPath());
					r.addSource(getClass().getResourceAsStream("/resources/netMHC.R"));
					r.run(true);
				} catch (Throwable e) {
					context.getLog().log(Level.SEVERE, "Could not plot!", e);
				}
			}
			return null;
		}
		
	}
	
	
	public static class PrepareFdrProgram extends GediProgram {

		public PrepareFdrProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.input);
			
			addOutput(params.fdrOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			LineWriter out = getOutputWriter(0);
			out.writeLine("Mode,Peptide length,Annotation,ALC,targets,decoys,ambiguous,Q1,Q2");
			out.close();
			
			return null;
		}
	}
	

	public static class ComputeFdrProgram extends GediProgram {

		private LocationAnnotation[] include;
		private String mode;

		public ComputeFdrProgram(FindGenomicPeptidesParameterSet params, GediParameter<File> in, GediParameter<File> out, String mode, LocationAnnotation... include) {
			this.include = include;
			this.mode = mode;
			
			addInput(in);
			addInput(params.fdrOut);
			addInput(params.plot);
			
			addInput(params.input);
			
			addOutput(out);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			File fdr = getParameter(1);
			boolean plot = getBooleanParameter(2);
			
			String prefix = getParameter(3);
			
			String input = csv.getPath();
			
			HashSet<String> allowed = EI.wrap(include).map(l->l.name()).set();
			
			context.getLog().info("Compute FDR cutoffs for "+mode);
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			
			int[] decoyCounter = {0,0}; // decoy only, both
			
			// [len,prio]->ALC,[target,decoy,both]
			HashMap<MutablePair<Integer,Integer>,int[][]> counter = new HashMap<>();

			
			for (String[][] af : iteratePeaksBlocks(input, header,null, l->allowed.contains(l[header.Item.get("Annotation")]),false).loop()) {
				int len = EI.seq(0, af.length).mapToInt(i->removeMod(af[i][header.Item.get("Peptide")]).length()).unique(false).getUniqueResult(true, true);
				int alc = EI.seq(0, af.length).mapToInt(i->Integer.parseInt(af[i][header.Item.get("ALC (%)")])).unique(false).getUniqueResult(true, true);

				int best = LocationAnnotation.valueOf(af[0][header.Item.get("Annotation")]).ordinal();
				for (int i=1; i<af.length; i++)
					best = Math.min(best,LocationAnnotation.valueOf(af[i][header.Item.get("Annotation")]).ordinal());

				int ubest = best;
				String[][] f = ArrayUtils.restrict(af, i->LocationAnnotation.valueOf(af[i][header.Item.get("Annotation")]).ordinal()==ubest);
		
				
				HashSet<String> decoy = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoy")]).set();
				if (decoy.contains("D"))
					decoyCounter[decoy.contains("T")?1:0]++;
				
				int ind = 0;
				if (decoy.contains("D"))
					ind=1;
				if (decoy.contains("D") && decoy.contains("T"))
					ind=2;
				
				counter.computeIfAbsent(new MutablePair<>(len,best), x->new int[3][101])[ind][alc]++;
				counter.computeIfAbsent(new MutablePair<>(len,-1), x->new int[3][101])[ind][alc]++;
			}
			
			context.getLog().info("Decoys unique: "+decoyCounter[0]);
			context.getLog().info("Decoy==target: "+decoyCounter[1]);
			double fac = 1/(decoyCounter[0]/((double)decoyCounter[0]+decoyCounter[1]));
			context.getLog().info("Decoy factor: "+fac);
			
			
			HashMap<MutablePair<Integer,Integer>,double[][]> qvals = new HashMap<>();
			for (MutablePair<Integer,Integer> l : counter.keySet()) {
				double[][] a;
				qvals.put(l, a = new double[2][101]);
				
				int[][] c = counter.get(l);
				ArrayUtils.cumSumInPlace(c[0], -1);
				ArrayUtils.cumSumInPlace(c[1], -1);
				ArrayUtils.cumSumInPlace(c[2], -1);
				
				for (int i=0; i<a[0].length; i++) {
					a[0][i] = c[1][i]/(double)c[0][i];
					if (Double.isNaN(a[0][i])) a[0][i] = 1;
					if (i>0) a[0][i] = Math.min(a[0][i], a[0][i-1]);
					
					a[1][i] = c[1][i]/((double)c[0][i]+(double)c[2][i])*fac;
					if (Double.isNaN(a[1][i])) a[1][i] = 1;
					if (i>0) a[1][i] = Math.min(a[1][i], a[1][i-1]);
				}
			}

			
			fdrLock.lock();
			context.getLog().info("Write FDR statistics for "+mode);
			LineWriter fout = new LineOrientedFile(fdr.getAbsolutePath()).append();
				
			for (MutablePair<Integer,Integer> pair : counter.keySet()) {
				int[][] m = counter.get(pair);
				double[][] a = qvals.get(pair);
				for (int i=0; i<m[0].length; i++) 
					fout.writef("%s,%d,%s,%d,%d,%d,%d,%3g,%3g\n", mode, pair.Item1, pair.Item2==-1?"-":LocationAnnotation.values()[pair.Item2], i, m[0][i], m[1][i], m[2][i], a[0][i], a[1][i]);
			}
			fout.close();
			
			
			fdrLock.unlock();
			
			context.getLog().info("Write peptide FDR list for "+mode);
			MutablePair<Integer,Integer> key = new MutablePair<>();
			LineWriter out = getOutputWriter(0);
			for (String[][] af : iteratePeaksBlocks(input, header,
													h->out.writeLine2(h+",UniqueAnnotation,Q1,Q2,aQ1,aQ2"), 
													l->allowed.contains(l[header.Item.get("Annotation")]),
													false
													).loop()) {
				
				int len = EI.seq(0, af.length).mapToInt(i->removeMod(af[i][header.Item.get("Peptide")]).length()).unique(false).getUniqueResult(true, true);
				int alc = EI.seq(0, af.length).mapToInt(i->Integer.parseInt(af[i][header.Item.get("ALC (%)")])).unique(false).getUniqueResult(true, true);

				int best = LocationAnnotation.valueOf(af[0][header.Item.get("Annotation")]).ordinal();
				for (int i=1; i<af.length; i++)
					best = Math.min(best,LocationAnnotation.valueOf(af[i][header.Item.get("Annotation")]).ordinal());

				int ubest = best;
				String[][] f = ArrayUtils.restrict(af, i->LocationAnnotation.valueOf(af[i][header.Item.get("Annotation")]).ordinal()==ubest);
				
				HashSet<String> decoy = EI.seq(0, f.length).map(i->f[i][header.Item.get("Decoy")]).set();
				
				String d = "T";
				if (decoy.contains("D"))
					d="D";
				if (decoy.contains("D") && decoy.contains("T"))
					d="B";
				
				key.Item1 = len;
				key.Item2 = best;
				double q1a = decoy.size()==1?qvals.get(key)[0][alc]:Double.NaN;
				double q2a = qvals.get(key)[1][alc];
				
				key.Item2 = -1;
				double q1 = decoy.size()==1?qvals.get(key)[0][alc]:Double.NaN;
				double q2 = qvals.get(key)[1][alc];
				

				String se = EI.seq(0, f.length).map(i->f[i][header.Item.get("Sequence")]).concat(";");
				String loc = EI.seq(0, f.length).map(i->f[i][header.Item.get("Location")]).concat(";");
				String ann = EI.seq(0, f.length).map(i->f[i][header.Item.get("Annotation")]).concat(";");
				String mann = EI.seq(0, f.length).map(i->LocationAnnotation.valueOf(f[i][header.Item.get("Annotation")])).sort().unique(true).concat(";");
				
				f[0][header.Item.get("Sequence")]=se;
				f[0][header.Item.get("Decoy")]=d;
				f[0][header.Item.get("Location")]=loc;
				f[0][header.Item.get("Annotation")]=ann;
				
				
				
//				for (int i=0; i<f.length; i++)
					out.writef("%s,%s,%.3g,%.3g,%.3g,%.3g\n",StringUtils.concat(",", f[0]),mann,q1,q2,q1a,q2a);
			}
			
			out.close();

			if (plot) {
				try {
					context.getLog().info("Running R script for plotting");
					RRunner r = new RRunner(prefix+".fdr"+mode+".R");
					r.set("file",getOutputFile(0).getPath());
					r.addSource(getClass().getResourceAsStream("/resources/fdr_types.R"));
					r.run(true);
				} catch (Throwable e) {
					context.getLog().log(Level.SEVERE, "Could not plot!", e);
				}
			}
			return null;
		}
		
	}
	
	
	public static class ToPeptideListProgram extends GediProgram {

		private LocationAnnotation[] include;

		public ToPeptideListProgram(FindGenomicPeptidesParameterSet params, GediParameter<File> out, LocationAnnotation... include) {
			this.include = include;
			
			addInput(params.annotatedPeaksOut);
			addInput(params.deltaNext);
			
			addOutput(out);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			File csv = getParameter(0);
			int minDeltaNext = getIntParameter(1);
			
			String input = csv.getPath();
			
			HashSet<String> allowed = EI.wrap(include).map(l->l.name()).set();
			
			context.getLog().info("Determine best PSM for "+EI.wrap(include).map(l->l.name()).concat(","));
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			HashMap<String,int[]> best = new HashMap<>();
			
			for (String[][] f : iteratePeaksBlocks(input, header,null, l->allowed.contains(l[header.Item.get("Annotation")]),false).loop()) {
				for (int i=0; i<f.length; i++) {
					String pep = f[i][header.Item.get("Sequence")];
					int alc = Integer.parseInt(f[i][header.Item.get("ALC (%)")]);
					int[] pp = best.computeIfAbsent(pep, x->new int[2]);
					pp[0] = Math.max(pp[0],alc);
					pp[1]++;
				}
			}
			
			
			
			context.getLog().info("Write peptide list for "+EI.wrap(include).map(l->l.name()).concat(","));
			LineWriter out = getOutputWriter(0);
			
			for (String[][] f : iteratePeaksBlocks(input, header,
													h->out.writeLine2(h+",Delta next,spectra"), 
													l->allowed.contains(l[header.Item.get("Annotation")]),
													false
													).loop()) {
				
				// identify block with same alc
				int e;
				for (e=1; e<f.length && Integer.parseInt(f[e-1][header.Item.get("ALC (%)")])==Integer.parseInt(f[e][header.Item.get("ALC (%)")]); e++);
				
				
				// this is not reasonable: if there is only the I/L difference between two sequences (i.e. both occurs in the target set of sequences), then why not keep both?
				boolean allSequencesSame = true;
				boolean anySequenceHasBetterSpectrum = false;
				int[] bestHit = null;
				for (int i=0; i<e; i++) {
					if (i>0 && !areEqualUptoIsobar(f[i-1][header.Item.get("Sequence")],f[i][header.Item.get("Sequence")]))
						allSequencesSame = false;
					bestHit = best.get(f[i][header.Item.get("Sequence")]);
					if (bestHit==null || bestHit[0]!=Integer.parseInt(f[i][header.Item.get("ALC (%)")]))
						anySequenceHasBetterSpectrum = true;
				}
				
				if (!allSequencesSame)
					continue;
				
				// remove it from the map (important if there is another spectrum with equal score)!
				if (!anySequenceHasBetterSpectrum) 
					best.remove(f[0][header.Item.get("Sequence")]);
				
				int deltaNext = e==f.length?100:(Integer.parseInt(f[0][header.Item.get("ALC (%)")])-Integer.parseInt(f[e][header.Item.get("ALC (%)")]));
				
				if (deltaNext>=minDeltaNext && !anySequenceHasBetterSpectrum) {
					for (int i=0; i<e; i++) {
						out.write(StringUtils.concat(",", f[i]));
						out.writef(",%d,%d\n",deltaNext,bestHit[1]);
					}
				}
				
			}
			
			out.close();

			return null;
		}
	}

	public static class FindGenomicPeptidesProgram extends GediProgram {



		public FindGenomicPeptidesProgram(FindGenomicPeptidesParameterSet params) {
			addInput(params.input);
			addInput(params.seq);
			addInput(params.nthreads);
			addInput(params.minlen);
			addInput(params.maxlen);
			addInput(params.deltaFirst);
			addInput(params.test);
			addInput(params.rnd);
			addInput(params.keepAll);
			addInput(params.spliced);
			addInput(params.deltaNext);
			
			addOutput(params.annotatedPeaksOut);
			addOutput(params.unidentifiedPeaksOut);
		}

		@Override
		public String execute(GediProgramContext context) throws Exception {
			String input = getParameter(0);
			String[] fasta = getParameters(1).toArray(new String[0]);
			int nthreads = getIntParameter(2);
			int minlen = getIntParameter(3);
			int maxlen = getIntParameter(4);
			int deltaFirst = getIntParameter(5);
			boolean test = getBooleanParameter(6);
			boolean rnd = getBooleanParameter(7);
			boolean all = getBooleanParameter(8);
			boolean spliced = getBooleanParameter(9);
			int deltaNext = getIntParameter(10);
			
			Genomic genomic;
			
			if (Genomic.check(fasta[0])!=null) {
				genomic = Genomic.get(fasta);
				fasta = Genomic.get(fasta).getGenomicFastaFiles().toArray(String.class);
				context.getLog().info("Recognized genomic input. Will annotate output!");
			} else {
				genomic = null;
				context.getLog().info("Recognized fasta. Will only search for peptides!");
			}
			
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = createTrie(context,input,minlen,maxlen);
			
			context.getLog().info("Preparing for Aho-Corasick");
			aho.prepareAhoCorasick();
			context.getLog().info("Trie has "+aho.getNodeCount()+" nodes!");
			
			
			if (test)
				searchTest(context, genomic, fasta, aho, nthreads, minlen,maxlen, rnd, spliced);
			else
				searchAll(context,genomic,fasta,aho,nthreads,minlen,maxlen, rnd, spliced);
			
			context.getLog().info("Writing output");
			
			LineWriter out = getOutputWriter(0);
			LineWriter out2 = getOutputWriter(1);
			ModifiedAminoAcid[] buff = new ModifiedAminoAcid[1024];
			for (int i=0; i<buff.length; i++) buff[i] = new ModifiedAminoAcid();
			int[] ind = new int[maxlen];
			char[] c = new char[maxlen];
			
			MutableMonad<HeaderLine> header = new MutableMonad<>(); 
			for (String[][] f : iteratePeaksBlocks(input, header, h->{
				out.writeLine2(h+",PSM rank,Location count,Delta first,Decoy,Location,Annotation,Sequence");
				out2.writeLine2(h+",Delta next");
			},null,true).loop()) {
				
				boolean found = false;
				int topalc = -1;
				for (int i=0; i<f.length; i++) {
					List<ImmutableReferenceGenomicRegion<String>> loc = new ArrayList<>();
					
					int alc = Integer.parseInt(f[i][header.Item.get("ALC (%)")]);
					int len = ModifiedAminoAcid.parse(f[i][header.Item.get("Peptide")],buff);
					
					if (len>=minlen && len<=maxlen) {
						for (int p=0; p<len; p++)
							c[p] = buff[p].aa;
						
						if (topalc==-1) topalc = alc;
						
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
							if (add!=null) 
								loc.addAll(add);
						}
						while (ArrayUtils.increment(ind, len,(p)->buff[p].getAlternative(buff,p)!=0?2:1));
					}
					
					
					if (loc!=null && loc.size()>0) {
						
						
						for (ImmutableReferenceGenomicRegion<String> l : loc) {
							if (topalc-alc<=deltaFirst || all) {
								out.write(StringUtils.concat(",", f[i]));
								out.writef(",%d,%d,%d,%s,%s,%s,%s\n",i+1,loc.size(),topalc-alc,isDecoy(l)?"D":"T",l.toLocationString(),annotateLocation(genomic, l).toString(),l.getData());
							}
						}
						found = true;
						
					} else if (all) {
						out.write(StringUtils.concat(",", f[i]));
						out.writef(",,,,,,,\n");
					}
				}
				
				if (!found) {
					int delta = f.length==1?100:Integer.parseInt(f[0][header.Item.get("ALC (%)")])-Integer.parseInt(f[1][header.Item.get("ALC (%)")]);
					String pep = removeMod(f[0][header.Item.get("Peptide")]);
					if (/*delta>=deltaNext && */pep.length()>=minlen && pep.length()<=maxlen) {
						out2.write(StringUtils.concat(",", f[0]));
						out2.writef(",%d\n",delta);
					}
				}
				
				
			}
			out.close();
			out2.close();

			return null;
		}
		

		private void searchAll(GediProgramContext context, Genomic genomic, String[] fasta,
				Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, int nthreads, int minlen, int maxlen, boolean rnd, boolean spliced) {
			
			FastaHeaderParser pars = new DefaultFastaHeaderParser(' ');
			
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
				.iff(genomic!=null && spliced, ei->
					genomic.getTranscripts().ei()
						.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Spliced peptides...")
						.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
						.map(tr->tr.getData().getCds(tr))
						.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
						.filter(fe->fe.getSequence().length()>=minlen*3)
						.unfold(fe->EI.wrap(new MutablePair<>(fe,42),new MutablePair<>(fe,43)))
						.chain(ei))
				.parallelized(nthreads, 1, ()->(rnd?new RandomNumbers():null), (ei,ran)->ei
				.map(ss->{
					String seq = ss.Item1.getSequence();
					String name = pars.getId(ss.Item1.getHeader());
					int frame = ss.Item2;
					if (isDna(seq))
						return searchDna(name,seq,frame,ran,aho,maxlen);
					else 
						return searchProtein(name,seq,frame==1,ran,aho);
				}))
				.iff(genomic==null, ei->ei.map(a->null))
				.removeNulls()
				.log(context.getLog(),Level.INFO);
			
			
			for (List<ImmutableReferenceGenomicRegion<String>> v : aho.values()) {
				HashSet<ImmutableReferenceGenomicRegion<String>> set = new HashSet<>(v);
				v.clear();
				v.addAll(set);
			}
			
		}
		
		private void searchTest(GediProgramContext context, Genomic genomic, String[] fasta,
				Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, int nthreads, int minlen, int maxlen, boolean rnd, boolean spliced) {
			
			FastaHeaderParser pars = new DefaultFastaHeaderParser(' ');
			
			context.getLog().info("Starting Aho-Corasick using "+nthreads+" threads");
			EI.singleton(new FastaEntry("22",genomic.getSequence(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain("22+"), new ArrayGenomicRegion(0,(int)genomic.getLength("22")))).toString()))
				.sideEffect(fe->{
					if (pars.getId(fe.getHeader()).contains(":") && genomic!=null)
						throw new RuntimeException("Sequence names may not contain : !");
				})
				.iff(genomic!=null, ei->
					genomic.getTranscripts().ei("22+").chain(genomic.getTranscripts().ei("22-"))
						.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
						.filter(fe->fe.getSequence().length()>=minlen*3)
						.chain(ei))
//				.head(1)
				.unfold(fe->{
					if (isDna(fe.getSequence()))
						return EI.seq(0, 12).map(frame->new MutablePair<>(fe, frame));
					else
						return EI.seq(0, 1).map(frame->new MutablePair<>(fe, frame));
				})
				.iff(genomic!=null && spliced, ei->
					genomic.getTranscripts().ei("22+")
						.progress(context.getProgress(), (int)genomic.getTranscripts().size(), x->"Spliced peptides...")
						.filter(tr->SequenceUtils.checkCompleteCodingTranscript(genomic, tr))
						.map(tr->tr.getData().getCds(tr))
						.map(tr->new FastaEntry(tr.toLocationString(), genomic.getSequence(tr).toString()))
						.filter(fe->fe.getSequence().length()>=minlen*3)
						.unfold(fe->EI.wrap(new MutablePair<>(fe,42),new MutablePair<>(fe,43)))
						.chain(ei))
				.parallelized(nthreads, 1, ()->(rnd?new RandomNumbers():null), (ei,ran)->ei
				.map(ss->{
					String seq = ss.Item1.getSequence();
					String name = pars.getId(ss.Item1.getHeader());
					int frame = ss.Item2;
					if (isDna(seq))
						return searchDna(name,seq,frame,ran,aho,maxlen);
					else 
						return searchProtein(name,seq,frame==0,ran,aho);
				}))
				.iff(genomic==null, ei->ei.map(a->null))
				.removeNulls()
				.log(context.getLog(),Level.INFO);
			
			
			for (List<ImmutableReferenceGenomicRegion<String>> v : aho.values()) {
				HashSet<ImmutableReferenceGenomicRegion<String>> set = new HashSet<>(v);
				v.clear();
				v.addAll(set);
			}
			
		}


		private boolean isDna(String seq) {
			return new Alphabet("ACGTN".toCharArray()).isValid(seq);
		}
		
		// frame 0-6: target, 6-12 decoy
		// frame odd: reverse
		private String searchDna(String name, String seq, int frame, RandomNumbers rnd,
				Trie<List<ImmutableReferenceGenomicRegion<String>>> aho, int maxl) {
		
			boolean transcr = name.contains(":");
			boolean decoy = frame>=6;
			boolean reverse = frame%2==1;
			ImmutableReferenceGenomicRegion<Object> parent = transcr?ImmutableReferenceGenomicRegion.parse(name):null;
			
			boolean spliced = frame>=42;
			if (spliced) {
				decoy=frame==43;
				CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, 0, false);
				if (decoy && rnd!=null)
					aa = rnd.shuffle(aa.toString());
				else if (decoy)
					aa = StringUtils.reverse(aa);

				findSpliced(aa.toString(),aho,parent,decoy,maxl);
				return null;
			}
			
			
			frame = (frame%6)/2;
			
			if (transcr && reverse) return null;
			
//			CharSequence aa = SequenceUtils.translate((!reverse?seq:SequenceUtils.getDnaReverseComplement(seq)).substring(frame));
			CharSequence aa = new SequenceUtils.SixFrameTranslatedSequence(seq, frame, reverse);
			if (decoy && rnd!=null)
				aa = rnd.shuffle(aa.toString());
			else if (decoy)
				aa = StringUtils.reverse(aa);
			
			int c = 0;
			for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(aa).loop()) {
				int s;
				if (!decoy)
					s = res.getStart()*3+frame;
				else
					s = (aa.length()-res.getEnd())*3+frame;
				
				if (reverse) s = seq.length()-s-res.getLength()*3;
				if (transcr) {
					GenomicRegion reg = parent.map(new ArrayGenomicRegion(s,s+res.getLength()*3));
					if (reg.getNumParts()>1) 
						res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()));
					// otherwise it's in the genome!
				} else {
					res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,!reverse), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()));
				}
				c++;
			}
			if (transcr)
				return null;
			return "For sequence "+name+" frame "+frame+" reverse "+reverse+" decoy "+decoy+": "+c;
		}
		
		private void findSpliced(String aa, Trie<List<ImmutableReferenceGenomicRegion<String>>> aho,
				ImmutableReferenceGenomicRegion<Object> parent, boolean decoy, int maxl) {
			
			int inter = 25;
			SplitSequence split = new SplitSequence(aa, 0, 1, 2, maxl-1);
			
			for (int s=1; s<aa.length()-2; s++) {
				split.len1 = Math.min(s, maxl-1);
				split.start1 = s-split.len1;
				for (int l=1; l<inter && s+l<aa.length()-1; l++) {
					split.start2=s+l;
					split.len2 = Math.min(maxl-1, aa.length()-s-l);
					
					
					for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(split).loop()) {
						if (res.getStart()<split.len1 & res.getEnd()>split.len1)  {
							
							GenomicRegion reg = new ArrayGenomicRegion(
										3*(res.getStart()+split.start1),
										3*(split.start1+split.len1),
										3*(split.start2),
										3*(split.start2+res.getLength()-split.len1+res.getStart())
										);
							if (decoy)
								reg = reg.reverse(parent.getRegion().getTotalLength());
							reg = parent.map(reg);
							
							if (reg.getNumParts()>1) 
								res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+parent.getReference().getName():parent.getReference().getName(),parent.getReference().getStrand()), reg, res.getKey().toString()));
							
						}
						
					}

				}
					
			}
			
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

		private String searchProtein(String name, String seq, boolean decoy, RandomNumbers rnd,
				Trie<List<ImmutableReferenceGenomicRegion<String>>> aho) {
		
			if (decoy && rnd!=null)
				seq = rnd.shuffle(seq);
			else if (decoy)
				seq = StringUtils.reverse(seq).toString();
			
			int c = 0;
			for (AhoCorasickResult<List<ImmutableReferenceGenomicRegion<String>>> res :  aho.iterateAhoCorasick(seq).loop()) {
				int s = res.getStart();
				if (decoy) s = seq.length()-s-res.getLength()*3;
				res.getValue().add(new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(decoy?"REV_"+name:name,true), new ArrayGenomicRegion(s,s+res.getLength()*3), res.getKey().toString()));
				c++;
			}
			return "For sequence "+name+(decoy?"decoy":"")+": "+c;
		}

		private boolean isDecoy(ImmutableReferenceGenomicRegion<String> hit) {
			return hit.getReference().getName().startsWith("REV");
		}
		
		private LocationAnnotation annotateLocation(Genomic g, ImmutableReferenceGenomicRegion<String> ahit) {
			if (g==null) return LocationAnnotation.CDS;
			
			if (ahit.getReference().getName().startsWith("REV_"))
				ahit = new ImmutableReferenceGenomicRegion<>(Chromosome.obtain(StringUtils.removeHeader(ahit.getReference().getName(), "REV_"),ahit.getReference().getStrand()), ahit.getRegion());
			
			ImmutableReferenceGenomicRegion<String> hit = ahit;
			ArrayList<ImmutableReferenceGenomicRegion<Transcript>> trs = g.getTranscripts().ei(hit).list();
			if (EI.wrap(trs)
				.filter(tr->tr.getData().isCoding())
				.map(tr->tr.getData().getCds(tr))
				.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()) && tr.getRegion().induce(hit.getRegion()).getStart()%3==0)
				.count()>0)
				return LocationAnnotation.CDS;
			
			if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding())
					.map(tr->tr.getData().get5Utr(tr))
					.filter(tr->tr.getRegion().intersects(hit.getRegion()) && tr.getRegion().isIntronConsistent(hit.getRegion()))
					.count()>0)
					return LocationAnnotation.UTR5;
			
			if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding())
					.map(tr->tr.getData().get3Utr(tr))
					.filter(tr->tr.getRegion().intersects(hit.getRegion()) && tr.getRegion().isIntronConsistent(hit.getRegion()))
					.count()>0)
					return LocationAnnotation.UTR3;
			
			
			if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding())
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
					.count()>0)
					return LocationAnnotation.OffFrame;
			
			if (EI.wrap(trs)
					.filter(tr->!tr.getData().isCoding())
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
					.count()>0)
					return LocationAnnotation.ncRNA;
			
			if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding())
					.map(tr->tr.getData().getCds(tr))
					.filter(tr->tr.getRegion().contains(hit.getRegion()))
					.count()>0)
					return LocationAnnotation.PeptideSpliced;
			
			
			if (EI.wrap(trs)
					.filter(tr->tr.getData().isCoding())
					.filter(tr->tr.getRegion().containsUnspliced(hit.getRegion()))
					.count()>0)
					throw new RuntimeException(hit.toLocationString());
			
			
			if (g.getGenes().ei(hit)
					.count()>0)
					return LocationAnnotation.Intronic;
			
			return LocationAnnotation.Intergenic;
		}
		
		private Trie<List<ImmutableReferenceGenomicRegion<String>>> createTrie(GediProgramContext context, String input, int minlen, int maxlen) throws IOException {
			context.getLog().info("Creating keyword trie");
			
			MutableInteger pepCol = new MutableInteger();
			ExtendedIterator<String> lit = EI.lines(input)
					.progress(context.getProgress(), -1, s->s)
					.skip(1,h->pepCol.N=new HeaderLine(h,',').get("Peptide"))
					.map(l->StringUtils.splitField(l, ',', pepCol.N));
			
			
			ModifiedAminoAcid[] buff = new ModifiedAminoAcid[1024];
			for (int i=0; i<buff.length; i++) buff[i] = new ModifiedAminoAcid();
			int[] ind = new int[maxlen];
			char[] c = new char[maxlen];
			
			Trie<List<ImmutableReferenceGenomicRegion<String>>> aho = new Trie<>();
			for (String l : lit.loop()) {
				
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
						List<ImmutableReferenceGenomicRegion<String>> list = Collections.synchronizedList(new ArrayList<ImmutableReferenceGenomicRegion<String>>());
						aho.put(kw, list);
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
			
			return aho;
		}
		
	}
	
	public static class ModifiedAminoAcid {
		char aa;
		String delta;
		
		public void set(char aa, String delta) {
			this.aa = aa;
			this.delta = delta;
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
		public GediParameter<String> seq = new GediParameter<String>(this,"seq", "DNA fasta", true, new StringParameterType());
		public GediParameter<String> tis = new GediParameter<String>(this,"tis", "Start codons to consider (prioritized)", true, new StringParameterType(),"AUG,CUG,ACG,GUG,AUC");

		public GediParameter<String> hla = new GediParameter<String>(this,"hla", "File containing HLA allels (each line); could also be in ${input%.csv}.hla, in which case you do not have to specify this parameter!", true, new StringParameterType(),true);
		public GediParameter<String> netmhc = new GediParameter<String>(this,"netmhc", "Command to call netmhc", true, new StringParameterType(),"netMHCpan");
		
		public GediParameter<Boolean> test = new GediParameter<Boolean>(this,"test", "Only search against chromosome 22", false, new BooleanParameterType());
		public GediParameter<Boolean> rnd = new GediParameter<Boolean>(this,"rnd", "Randomize decoys", false, new BooleanParameterType());

		public GediParameter<Boolean> keepAll = new GediParameter<Boolean>(this,"all", "Keep all lines", false, new BooleanParameterType());
		public GediParameter<Boolean> spliced = new GediParameter<Boolean>(this,"spliced", "Also search for spliced peptides", false, new BooleanParameterType());

		public GediParameter<File> annotatedPeaksOut = new GediParameter<File>(this,"${in}.annotated.csv.gz", "Peaks file annotated with matched sequence, location and type (filtered by df and pep length)", false, new FileParameterType());
		public GediParameter<File> unidentifiedPeaksOut = new GediParameter<File>(this,"${in}.unidentified.csv.gz", "Peaks file of all top unidentified hits with delta next filter applied", false, new FileParameterType());
		
		public GediParameter<File> unidentifiedPeaksFOut = new GediParameter<File>(this,"${in}.pep.unidentified.csv", "Peptide list of unannotated spectra", false, new FileParameterType());
		
		public GediParameter<File> pepCdsOut = new GediParameter<File>(this,"${in}.pep.CDS.tmp", "Peaks file filtered for CDS hits and best hit per Sequence (filtered by dn)", false, new FileParameterType()).setRemoveFile(true);
		public GediParameter<File> pepTranscriptomeOut = new GediParameter<File>(this,"${in}.pep.transcriptome.tmp", "Peaks file filtered for CDS/UTR/OffFrame/ncRNA hits and best hit per Sequence (filtered by dn)", false, new FileParameterType()).setRemoveFile(true);
		public GediParameter<File> pepGenomeOut = new GediParameter<File>(this,"${in}.pep.genome.tmp", "Peaks file filtered for best hit per Sequence (filtered by dn)", false, new FileParameterType()).setRemoveFile(true);

		public GediParameter<File> fdrOut = new GediParameter<File>(this,"${in}.fdr.csv", "FDR statistics", false, new FileParameterType());
		public GediParameter<File> fdrPlot = new GediParameter<File>(this,"${in}.fdr.pdf", "FDR statistics plot", false, new FileParameterType());

		public GediParameter<File> pepCdsFOut = new GediParameter<File>(this,"${in}.pep.CDS.csv", "Peaks file filtered for CDS hits and best hit per Sequence (filtered by dn)", false, new FileParameterType());
		public GediParameter<File> pepTranscriptomeFOut = new GediParameter<File>(this,"${in}.pep.transcriptome.csv", "Peaks file filtered for CDS/UTR/OffFrame/ncRNA hits and best hit per Sequence (filtered by dn)", false, new FileParameterType());
		public GediParameter<File> pepGenomeFOut = new GediParameter<File>(this,"${in}.pep.genome.csv", "Peaks file filtered for best hit per Sequence (filtered by dn)", false, new FileParameterType());
		
		public GediParameter<File> pepCdsAOut = new GediParameter<File>(this,"${in}.pep.CDS.annotated.csv", "Annotated peaks file (Orfs, NetMHC binding)", false, new FileParameterType());
		public GediParameter<File> pepTranscriptomeAOut = new GediParameter<File>(this,"${in}.pep.transcriptome.annotated.csv", "Annotated peaks file (Orfs, NetMHC binding)", false, new FileParameterType());
		public GediParameter<File> pepGenomeAOut = new GediParameter<File>(this,"${in}.pep.genome.annotated.csv", "Annotated peaks file (Orfs, NetMHC binding)", false, new FileParameterType());

		
		public GediParameter<Boolean> plot = new GediParameter<Boolean>(this,"plot", "Use R to produce various plots",false, new BooleanParameterType());

}


	
	
	
}
