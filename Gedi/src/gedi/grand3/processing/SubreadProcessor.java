package gedi.grand3.processing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.function.Supplier;
import java.util.function.UnaryOperator;
import java.util.logging.Logger;

import gedi.centeredDiskIntervalTree.CenteredDiskIntervalTreeStorage;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsVariation;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.reference.ReferenceSequence;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.knmatrix.SubreadCounterKNMatrices;
import gedi.grand3.knmatrix.SubreadCounterKNMatrixPerTarget;
import gedi.grand3.reads.ReadSource;
import gedi.grand3.targets.SnpData;
import gedi.grand3.targets.TargetCollection;
import gedi.util.SequenceUtils;
import gedi.util.datastructure.collections.intcollections.IntArrayList;
import gedi.util.datastructure.mapping.OneToManyMapping;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.CompoundParallelizedState;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;
import gedi.util.userInteraction.progress.Progress;

public class SubreadProcessor<A extends AlignedReadsData>  {

	private int nthreads = 0;
	private Genomic genomic;
	private SnpData masked;
	
	private boolean debug = false;
	
	private ReadSource<A> source;
	private Logger logger;

	
	
	public SubreadProcessor(Genomic genomic, ReadSource<A> source, 
			SnpData masked, Logger logger) {
		this.genomic = genomic;
		this.source = source;
		this.masked = masked;
		this.logger = logger;
	}
	
	public SubreadProcessor<A>  setNthreads(int nthreads) {
		this.nthreads = nthreads;
		return this;
	}
	
	public SubreadProcessor<A> setDebug(boolean debug) {
		this.debug = debug;
		return this;
	}
	
	
	public void writeSubreads(Supplier<Progress> progress, 
			TargetCollection targets, File output, ExperimentalDesign design, String[] subreads, Resimulator resimulator, String[] newConditions, int[][] conditionMapping) throws IOException {
		
		ExtendedIterator<ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>> srit = targets.iterateRegions()
				//EI.wrap(targets.getRegion("MYC"))
			.iff(progress!=null,ei->ei.progress(progress.get(), targets.getNumRegions(), r->"Processing "+r.getData()))
			.parallelizedState(nthreads, 5, resimulator==null?null:resimulator.createState(),(ei,resim)->ei.unfold(target->{
				
				ImmutableReferenceGenomicRegion<String>  currentTargetExtended = new ImmutableReferenceGenomicRegion<>(
						target.getReference(), 
						target.getRegion().extendAll(100000, 100000).intersect(0, genomic.getLength(target.getReference().getName())),
						target.getData());
				char[] targetSeq = genomic.getSequence(currentTargetExtended).toString().toUpperCase().toCharArray();
				resim.setSequence(currentTargetExtended,targetSeq);
				
				return source.getSubReads(target,null).map(r->{
					targets.classify(target, r, source.getStrandness(), true, resim);
					return resimulator.resimulate(resim,r);
				}).removeNulls();
			}))
			.iff(progress!=null,ei->ei.progress(progress.get(), -1, r->"Finished subreads"));
		

		if (newConditions!=null) {
			srit = srit.map(r->{
				SubreadsAlignedReadsData rd = r.getData().selectMergeConditions(newConditions.length, conditionMapping );
				if (rd==null) return null;
				return new ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData>(r.getReference(),r.getRegion(),rd);
			}).removeNulls();		
		}
		
		CenteredDiskIntervalTreeStorage<SubreadsAlignedReadsData> scit = new CenteredDiskIntervalTreeStorage<>(output.getPath(),SubreadsAlignedReadsData.class);
		scit.fill(srit,progress.get());
		
		DynamicObject cond = newConditions!=null?DynamicObject.from("conditions", DynamicObject.arrayOfObjects("name", newConditions)):DynamicObject.from("conditions", DynamicObject.arrayOfObjects("name", EI.seq(0, design.getCount()).map(index->design.getFullName(index))));
		DynamicObject subr = DynamicObject.from("subreads", DynamicObject.from(subreads));
		scit.setMetaData(DynamicObject.merge(cond,subr));

	}


	public void process(Supplier<Progress> progress, 
			TargetCollection targets,
			SubreadCounter... counter) {
		
		process(progress,
				targets,new OneToManyMapping<String,String>(),
				new DummyTargetCounter(),
				counter)
		.drain();
		
	}
	
	public <T extends TargetCounter<T,R>,R> ExtendedIterator<R> process(Supplier<Progress> progress, 
			TargetCollection targets, OneToManyMapping<String,String> mapper,
			T res,
			SubreadCounter... counter) {
		
		CompoundParallelizedState state = new CompoundParallelizedState();
		state.add(new SubreadProcessorMismatchBuffer(source.getConverter().getSemantic().length));
		state.add((SubreadCounter)res);
		for (SubreadCounter c : counter)
			state.add(c);
		
		ExtendedIterator<? extends List<ImmutableReferenceGenomicRegion<String>>> it;
		int num;
		if (res.mergeWithSameName()) {
			logger.info("Checking target names...");
			HashMap<String, ArrayList<ImmutableReferenceGenomicRegion<String>>> merge = targets.iterateRegions().indexMulti(r->mapper.inverseOrDefault(r.getData(),r.getData()));
			num = merge.size();
			it = EI.wrap(merge.values());
			
			if (num!=targets.getNumRegions())
				logger.info(String.format("...from %d targets, %d remain after merging by name!",targets.getNumRegions(),num));
			else
				logger.info(String.format("...all %d names are unique!",targets.getNumRegions()));

		} else {
			it = targets.iterateRegions().map(r->Arrays.asList(r));
			num = targets.getNumRegions();
		}
		
//		it = EI.singleton(Arrays.asList(genomic.getGeneMapping().apply("ENSG00000203326").toImmutable()));		
		return it
			.iff(progress!=null,ei->ei.progress(progress.get(), num, r->"Processing "+r.get(0).getData()))
			.parallelizedState(nthreads, 5, state, (ei,b)->ei.map(targetset->{
				TargetCounter<T, R> t = b.get(1);
				t.startChunk();
				for (ImmutableReferenceGenomicRegion<String> target : targetset) {
					processTarget(targets,target,b);
				}
				return t.getResultsForCurrentTargets();
			}))
			.iff(progress!=null,ei->ei.progress(progress.get(), targets.getNumRegions(), r->"Finished"))
			.unfold(l->EI.wrap(l));
	}
	
	
	public void processTarget(Supplier<Progress> progress, String name, 
			TargetCollection targets,
			SubreadCounter res) throws IOException {
		
		CompoundParallelizedState state = new CompoundParallelizedState();
		state.add(new SubreadProcessorMismatchBuffer(source.getConverter().getSemantic().length));
		state.add(res);
		
		ImmutableReferenceGenomicRegion<String> target = targets.getRegion(name);
		if (target==null) throw new RuntimeException("Target with name "+name+" unknown!");
		processTarget(targets,target,state);
	}

	
	private <T extends SubreadCounter<T>> void processTarget(
			TargetCollection targets, 
			ImmutableReferenceGenomicRegion<String> target, 
			CompoundParallelizedState state) {
		
		SubreadProcessorMismatchBuffer buff = state.<SubreadProcessorMismatchBuffer>get(0);
		
		ImmutableReferenceGenomicRegion<String>  currentTargetExtended = new ImmutableReferenceGenomicRegion<>(
				target.getReference(), 
				target.getRegion().extendAll(100000, 100000).intersect(0, genomic.getLength(target.getReference().getName())),
				target.getData());
		char[] targetSeq = genomic.getSequence(currentTargetExtended).toString().toUpperCase().toCharArray();

		ReferenceSequence refInd = target.getReference().toStrandIndependent();
		
		for (ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> sread : source.getSubReads(target,null).loop()) {
			// classify read
			targets.classify(target, sread, source.getStrandness(), true, buff);
			
			if (buff.getMode().equals(ReadCountMode.No)) 
				continue;

			char[] readSeq;
			// obtain sequence
			if (currentTargetExtended.getRegion().contains(sread.getRegion())) {
				readSeq = SequenceUtils.extractSequence(currentTargetExtended.induce(sread.getRegion()), targetSeq);
				if (!sread.getReference().getStrand().equals(currentTargetExtended.getReference().getStrand()))
					SequenceUtils.getDnaReverseComplementInplace(readSeq);
			}
			else {
				readSeq = genomic.getSequence(sread).toString().toUpperCase().toCharArray();
			}
			
			int len = sread.getRegion().getTotalLength();
			
			SubreadsAlignedReadsData data = sread.getData();
			// count mismatches per distinct
			for (int d=0; d<data.getDistinctSequences(); d++) {
				buff.distinct = d;
				
				// count totals
				for (int s=0; s<data.getNumSubreads(d); s++) {
					for (int i=data.getSubreadStart(d, s); i<data.getSubreadEnd(d, s, len); i++) {
						int gpos = sread.map(i);
						boolean indeletion = sread.getData().mapToRead1(d, i)==-1;
						if (!indeletion && !masked.isSnp(refInd, gpos)) {
							buff.increment(data.getSubreadId(d, s),readSeq[i]);
						}
					}
				}
				
				for (int v=0; v<data.getVariationCount(d); v++) {
					if (data.isMismatch(d, v)) {
						int gpos = sread.map(data.getMismatchPos(d, v));
						if (!masked.isSnp(refInd, gpos) && checkSequence(data.getMismatchGenomic(d, v).charAt(0),readSeq[data.getMismatchPos(d, v)],sread,readSeq,sread.getData().getVariation(d, v))) {
							int s = data.getSubreadIndexForPosition(d, data.getMismatchPos(d, v), len);
							buff.increment(data.getSubreadId(d, s), data.getMismatchGenomic(d, v).charAt(0), data.getMismatchRead(d, v).charAt(0));
						}
					}
				}
				
				if (debug) {
					System.out.println("Process:");
					System.out.println(d+" "+sread);
					System.out.println(buff);
					System.out.println();
				}
				
				for (int i=1; i<state.size(); i++)
					state.<SubreadCounter>get(i).count(buff);
				
				buff.reset();
			}
			
				
			
		}
		
	}

	
	
	private boolean checkSequence(char mm, char genome, ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> sread, char[] readSeq, AlignedReadsVariation var) {
		if (mm=='N' && genome!='A'&&genome!='C'&&genome!='G'&&genome!='T')
			return false;
		if (mm!=genome)
			throw new RuntimeException("Sequences do not match! This is a sign that references for read mapping and Grand3 are different!\n"+sread+"\nReference sequence: "+String.valueOf(readSeq)+"\n"+var);
		return true;

	}



	private static class DummyTargetCounter implements TargetCounter<DummyTargetCounter, Integer> {

		@Override
		public DummyTargetCounter spawn(int index) {
			return new DummyTargetCounter();
		}

		@Override
		public void integrate(DummyTargetCounter other) {
		}

		@Override
		public boolean mergeWithSameName() {
			return false;
		}
		
		@Override
		public void startChunk() {
		}

		@Override
		public void count(SubreadProcessorMismatchBuffer buffer) {
		}

		@Override
		public List<Integer> getResultsForCurrentTargets() {
			return Arrays.asList(1);
		}
		
		
	}
	
}
