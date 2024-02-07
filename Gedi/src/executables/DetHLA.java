package executables;

import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import cern.colt.bitvector.BitVector;
import gedi.proteomics.digest.FullAfterAADigester;
import gedi.util.StringUtils;
import gedi.util.datastructure.tree.Trie;
import gedi.util.functions.EI;
import gedi.util.io.text.HeaderLine;
import gedi.util.io.text.LineIterator;
import gedi.util.io.text.LineWriter;
import gedi.util.io.text.fasta.FastaEntry;
import gedi.util.io.text.fasta.FastaFile;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;


public class DetHLA {
	
	




	public static void main(String[] args) throws IOException {
		
		DetHLAParameterSet params = new DetHLAParameterSet();
		GediProgram pipeline = GediProgram.create("DetHLA",
				new SimulatePeptidesProgram(params),
				new InferHLAProgram(params)
				);
		GediProgram.run(pipeline, null, null, new CommandLineHandler("DetHLA","DetHLA.",args));

	}
	
	
	public static class InferHLAProgram extends GediProgram {

		
		public InferHLAProgram(DetHLAParameterSet params) {
			addInput(params.peptides);
			addInput(params.prefix);
			
			addOutput(params.out);
		}

		
		@Override
		public String execute(GediProgramContext context) throws Exception {
			File peps = getParameter(0);
			
			HashSet<String> peptides = new HashSet<>();
			
			LineIterator lit = EI.lines(peps.getPath());
			if (!lit.hasNext()) throw new RuntimeException("Peptide list is empty!");
			String first = lit.next();
			HeaderLine header = new HeaderLine(first,';');
			if (header.hasField("Peptide"))
				lit.splitField(';',header.get("Peptide")).map(s->removeMod(s)).toCollection(peptides);
			else
				EI.singleton(first).chain(lit).map(s->removeMod(s)).toCollection(peptides);

			
			Trie<ArrayList<String>> aho = new Trie<>();
			for (String p : peptides)
				aho.put(p, new ArrayList<>());
			
			HashMap<String, Integer> index = EI.wrap(peptides).indexPosition();
			
			HashMap<String,String> seqs = new FastaFile.EntryIterator(new InputStreamReader(getClass().getResourceAsStream("/resources/hla_prot.fasta")),false)
					.index(e->e.getHeader(),e->e.getSequence());
			HashMap<String,BitVector> header2NotPeps = new HashMap<>();
			
			HashMap<String,HashSet<String>> allele2Headers = new HashMap<>();
			allele2Headers.put("A", new HashSet<>());
			allele2Headers.put("B", new HashSet<>());
			allele2Headers.put("C", new HashSet<>());
			
			for (String h : seqs.keySet()) {
				
				String allele = getAllele(h).substring(0, 1);
				if (allele2Headers.containsKey(allele)) {
					
					BitVector found = new BitVector(index.size());
					aho.iterateAhoCorasick(seqs.get(h)).forEachRemaining(hit->{
						found.putQuick(index.get(hit.getKey().toString()),true);
						hit.getValue().add(h);
						});
					header2NotPeps.computeIfAbsent(getAllele(h,2),x->new BitVector(index.size())).or(found);
					allele2Headers.get(allele).add(getAllele(h,2));
				}
			}

			for (BitVector b : header2NotPeps.values())
				b.not();
			
			for (String a : allele2Headers.keySet()) {
				ArrayList<String> l = new ArrayList<>(allele2Headers.get(a));
				Collections.sort(l);
				
				ArrayList<int[]> best = new ArrayList<>();
				int binter = peptides.size();
				
				BitVector[] bvs = EI.wrap(l).map(h->header2NotPeps.get(h)).toArray(BitVector.class);
				for (int i=0; i<l.size(); i++) {
					BitVector buff = new BitVector(index.size());
					for (int j=i+1; j<l.size(); j++) {
						buff.replaceFromToWith(0, index.size()-1, bvs[i], 0);
						buff.and(bvs[j]);
						int le = buff.cardinality();
						if (le<=binter) {
							if (le<binter) best.clear();
							best.add(new int[] {i,j});
							binter = le;
						}
					}
				}
				System.out.println(a+" "+best.size());
				System.out.println(a+" "+binter+" "+EI.wrap(best).map(arr->l.get(arr[0])+","+l.get(arr[1])).concat(" "));
			}
			
			return null;
		}
		
		private static String removeMod(String s) {

			StringBuilder sb = new StringBuilder();
			int len = 0;
			int p = 0;
			for (int sepIndex=s.indexOf('('); sepIndex>=0; sepIndex = s.indexOf('(',p)) {
				for (int i=p; i<sepIndex-1; i++)
					sb.append(s.charAt(i));
				sb.append(s.charAt(sepIndex-1));
				
				p = s.indexOf(')',p)+1;
			}
			for (int i=p; i<s.length(); i++)
				sb.append(s.charAt(i));
			return sb.toString();
		}
	}

	
	
	public static class SimulatePeptidesProgram extends GediProgram {

		
		public SimulatePeptidesProgram(DetHLAParameterSet params) {
			addInput(params.n);
			addInput(params.hla);
			
			addInput(params.prefix);
			
			addOutput(params.peptides);
		}

		
		@Override
		public String execute(GediProgramContext context) throws Exception {
			
			int n = getIntParameter(0);
			File hla = getParameter(1);
			
			RandomNumbers rnd = new RandomNumbers();
			

			Trie<FastaEntry> full = new Trie<>();
			new FastaFile.EntryIterator(new InputStreamReader(getClass().getResourceAsStream("/resources/hla_prot.fasta")),false)
					.forEachRemaining(e->{
						full.put(getAllele(e.getHeader()), e);
					});
			
			ArrayList<String> hlas;
			if (hla!=null && hla.exists()) {
				hlas = EI.lines(hla.getPath()).map(s->StringUtils.removeHeader(s, "HLA-").replaceAll("\\*", "")).list();
			}
			else {
				ArrayList<String> all = new ArrayList<>(full.keySet());
				rnd.shuffle(all);
				hlas = new ArrayList<>();
				EI.wrap(all).filter(s->s.startsWith("A")).head(2).toCollection(hlas);
				EI.wrap(all).filter(s->s.startsWith("B")).head(2).toCollection(hlas);
				EI.wrap(all).filter(s->s.startsWith("C")).head(2).toCollection(hlas);
			}
			
			if (hlas.size()!=6) throw new RuntimeException("HLAs != 6!");
			if (EI.wrap(hlas).filter(h->full.getKeysByPrefix(h).size()>0).removeNulls().count()!=6) throw new RuntimeException("Sequences != 6!");

			
			FullAfterAADigester digest = new FullAfterAADigester(2, 'K','R');
			
			ArrayList<String> allpeps = EI.wrap(hlas)
				.map(h->{
					ArrayList<FastaEntry> l = new ArrayList<>(full.getValuesByPrefix(h));
					rnd.shuffle(l);
					return l.get(0);
				})
				.unfold(e->digest.iteratePeptides(e.getSequence()))
				.filter(p->p.length()>=8 && p.length()<=20)
				.unique(false)
				.list();
			rnd.shuffle(allpeps);
			
			LineWriter out = getOutputWriter(0);
			for (int i=0; i<n; i++)
				out.writeLine(allpeps.get(i));
			out.close();
			
			
			return null;
		}


		
	}
	
	public static String getAllele(String header) {
		return StringUtils.removeHeader(StringUtils.splitField(header, " ", 1), "HLA-").replaceAll("\\*", "");
	}
	
	public static String getAllele(String header, int acc) {
		String re = getAllele(header);
		for (int i=0; i<re.length(); i++) {
			if (re.charAt(i)==':') {
				if (acc--==1)
					return re.substring(0,i);
			}
		}

		return re;
	}
	
	
	
	public static class DetHLAParameterSet extends GediParameterSet {

		public GediParameter<String> prefix = new GediParameter<String>(this,"prefix", "Prefix for output files", true, new StringParameterType());

		public GediParameter<Integer> n = new GediParameter<Integer>(this,"npep", "Number of peptides to generate", false, new IntParameterType(),100);
		public GediParameter<File> hla = new GediParameter<File>(this,"${prefix}.hla", "", false, new FileParameterType(),true);
		
		
		public GediParameter<File> peptides = new GediParameter<File>(this,"${prefix}.csv.gz", "Peaks file or list of peptides", false, new FileParameterType());
		public GediParameter<File> out = new GediParameter<File>(this,"${prefix}.inferred", "Result", false, new FileParameterType());
		

	}

	
}
