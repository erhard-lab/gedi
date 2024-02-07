package gedi.grand10x.javapipeline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.logging.Level;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.BarcodedAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand10x.TrimmedGenomicRegion;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.tree.Trie;
import gedi.util.functions.EI;
import gedi.util.math.stat.descriptive.WeightedMeanVarianceOnline;
import gedi.util.mutable.MutableInteger;
import gedi.util.mutable.MutablePair;
import gedi.util.mutable.MutableTriple;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import gedi.util.r.RRunner;
import gedi.util.sequence.DnaSequence;

public class Grand10xClusterPositionInferenceProgram extends GediProgram {

	public Grand10xClusterPositionInferenceProgram(Grand10xParameterSet params) {
		addInput(params.prefix);
		addInput(params.reads);
		addInput(params.genomic);
		addInput(params.nthreads);
		addInput(params.plot);
		
		addOutput(params.infer3pparam);
		addOutput(params.infer3ptsv);
	}

	@Override
	public String execute(GediProgramContext context) throws Exception {
		String prefix = getParameter(0);
		GenomicRegionStorage<BarcodedAlignedReadsData> reads = getParameter(1);
		Genomic genomic = getParameter(2);
		int nthreads = getIntParameter(3);
		boolean plot = getParameter(4);

		
		String[] signals = {"AATAAA","ATTAAA","AGTAAA","TATAAA"};
		Trie<Integer> strie = new Trie<>();
		for (int i=0; i<signals.length; i++)
			strie.put(signals[i], i);
		
		int[] signalrange = {12,30}; // most of the ccds transcripts have the start of the aataaa in there! 
		int minDistOther = 1000;
		
		int[] count = new int[6];
		
		ArrayList<ImmutableReferenceGenomicRegion<Transcript>> niceTranscripts = genomic.getTranscripts().ei()
			.progress(context.getProgress(), (int)genomic.getTranscripts().size(), t->t.getData().getTranscriptId())
			.sideEffect(x->count[0]++)
			.filter(t->t.getData().isCoding())
			.sideEffect(x->count[1]++)
			.filter(t->{
				int t3p = GenomicRegionPosition.ThreePrime.position(t);
				int t3p1 = GenomicRegionPosition.ThreePrime.position(t,1);
				return genomic.getTranscripts()
					.ei(new ImmutableReferenceGenomicRegion<>(t.getReference(), new ArrayGenomicRegion(Math.min(t3p1,t3p),Math.max(t3p1,t3p))))
					.filter(o->!o.equals(t))
					.filter(o->Math.abs(GenomicRegionPosition.ThreePrime.position(o)-t3p)<minDistOther)
					.count()==0;
			})
			.sideEffect(x->count[2]++)
			.map(t->new MutableTriple<>(t,genomic.getSequence(t.getData().get3Utr(t)).toString(),genomic.getSequence(t.getDownstream(50)).toString()))
			.filter(tr->{
				String expectPas = tr.Item2.substring(Math.max(0, tr.Item2.length()-signalrange[1]));
				int paspos = expectPas.indexOf(signals[0]);
				return paspos>=0 && paspos<signalrange[1]-signalrange[0];
			})
			.sideEffect(x->count[3]++)
			.filter(tr->tr.Item2.length()>signalrange[1]?
				strie.iterateAhoCorasick(tr.Item2.substring(0,tr.Item2.length()-signalrange[1])).count()==0
				:true
			)
			.sideEffect(x->count[4]++)
			.filter(tr->
				SequenceUtils.getPolyAStretches(tr.Item2+tr.Item3).size()==0
			)
			.sideEffect(x->count[5]++)
			.map(tr->tr.Item1)
			.list();
		
		context.logf("Initial transcripts:  %d",count[0]);
		context.logf("Removed non-coding:   %d",count[1]);
		context.logf("Removed unclear 3':   %d",count[2]);
		context.logf("Removed no PAS:       %d",count[3]);
		context.logf("Removed another PAS:  %d",count[4]);
		context.logf("Removed pA stretches: %d",count[5]);
		
		Collections.shuffle(niceTranscripts);
		
		ClusterPositionStatistics stat = EI.wrap(niceTranscripts)
			.parallelized(nthreads, 8, ei->ei.map(t->{
				TrimmedGenomicRegion tgr = new TrimmedGenomicRegion();
				BarcodeCounter bco = new BarcodeCounter();
				for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : reads.ei(t).filter(r->tgr.set(r.getRegion()).isCompatibleWith(t.getRegion())).loop()) {
					for (int c=0; c<read.getData().getNumConditions(); c++)
						for (int d=0; d<read.getData().getDistinctSequences(); d++)
							bco.addBarcodes(read.getData().getBarcodes(d, c));
				}
				HashMap<DnaSequence, MutableInteger> bcoo = bco.get();
				ClusterPositionStatistics tmvo = new ClusterPositionStatistics();
				for (ImmutableReferenceGenomicRegion<BarcodedAlignedReadsData> read : reads.ei(t).filter(r->tgr.set(r.getRegion()).isCompatibleWith(t.getRegion())).loop()) {
					tgr.set(read.getRegion());
					int s = read.getReference().isMinus()?tgr.getStop():tgr.getStart();
					if (t.getRegion().contains(s)) {
						int dist = t.getRegion().getTotalLength()-t.induce(s);
						if (dist<=450) {
							double w = 0;
							for (int c=0; c<read.getData().getNumConditions(); c++)
								for (int d=0; d<read.getData().getDistinctSequences(); d++)
									for (DnaSequence bc : read.getData().getBarcodes(d, c)) 
										w+=1.0/bcoo.get(bc).N;
							tmvo.add(dist, w);
						}
					}
				}
				return new MutablePair<>(tmvo,t);
			}))
			.progress(context.getProgress(),(int)niceTranscripts.size(), t->t.Item2.toLocationString())
			.reduce(new ClusterPositionStatistics(),(a,b)->b.add(a.Item1));

		FileUtils.writeAllLines(new String[] {"Name\tValue","Mean\t"+stat.mvo.getMean(),"Sd\t"+stat.mvo.getStandardDeviation(),""}, getOutputFile(0));

		EI.seq(0, stat.histo.length).filterInt(i->stat.histo[i]>0).map(i->i+"\t"+stat.histo[i]).print("Distance\tFrequency", getOutputFile(1).getPath());
		
		
		if (plot) {
			try {
				context.getLog().info("Running R scripts for plotting");
				RRunner r = new RRunner(prefix+".plotclusterdist.R");
				r.set("prefix",prefix);
				r.addSource(getClass().getResourceAsStream("/resources/R/plotclusterdist.R"));
				r.run(true);
			} catch (Throwable e) {
				context.getLog().log(Level.SEVERE, "Could not plot!", e);
			}
		}
		
		return null;
	}
	
	
	private static class ClusterPositionStatistics {
		private WeightedMeanVarianceOnline mvo = new WeightedMeanVarianceOnline();
		private double[] histo = new double[1000];
		public void add(int d, double w) {
			mvo.add(d,w);
			histo[d]+=w;
		}
		
		
		public ClusterPositionStatistics add(ClusterPositionStatistics other) {
			this.mvo.add(other.mvo);
			ArrayUtils.add(this.histo, other.histo);
			return this;
		}
	}
	
	private static class BarcodeCounter {
		private HashMap<DnaSequence,MutableInteger> counter = new HashMap<DnaSequence,MutableInteger>();
		
		public void addBarcodes(DnaSequence[] dna) {
			for (DnaSequence d : dna) {
				counter.computeIfAbsent(d, x->new MutableInteger()).N++;
			}
		}
		
		public HashMap<DnaSequence, MutableInteger> get() {
			return counter;
		}
	}
}

