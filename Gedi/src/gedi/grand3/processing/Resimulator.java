package gedi.grand3.processing;

import java.util.HashMap;

import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.data.reads.SubreadsAlignedReadsData;
import gedi.core.genomic.Genomic;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.grand3.estimation.MismatchMatrix;
import gedi.grand3.estimation.ModelStructure;
import gedi.grand3.estimation.TargetEstimationResult.ModelType;
import gedi.grand3.estimation.models.Grand3Model;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.util.SequenceUtils;
import gedi.util.functions.EI;

public class Resimulator  {

	public enum ResimulatorModelType {
		None,Binom,TbBinom
	}
	

	private ResimulatorModelType modelType;
	private ModelStructure[][][] models;
	private HashMap<String, Double>[][][] ntrs;
	private MetabolicLabelType[] types;
	private MismatchMatrix mmMat;
	private long seed;
	private Genomic genomic;

	public Resimulator(Genomic genomic,ResimulatorModelType modelType, ModelStructure[][][] models, HashMap<String, Double>[][][] ntrs, MismatchMatrix mmMat, MetabolicLabelType[] types, long seed) {
		this.genomic = genomic;
		this.modelType = modelType;
		this.models = models;
		this.ntrs = ntrs;
		this.types = types;
		this.mmMat = mmMat;
		this.seed = seed;
		
		if (EI.wrap(types).map(t->t.getGenomic()).set().size()!=types.length)
			throw new RuntimeException("Simulation with multiple labling types only possible, when genomic nt are distinct!");
	}


	public ResimulateState createState() {
		return new ResimulateState(seed);
	}


	public ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> resimulate(ResimulateState state,
			ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> t) {
		
		if (state.getTargets()==null) return null;
		if (modelType.equals(ResimulatorModelType.None)) return t;
		
		char[] readSeq;
		// obtain sequence
		if (state.getSequenceRegion().getRegion().contains(t.getRegion())) {
			readSeq = SequenceUtils.extractSequence(state.getSequenceRegion().induce(t.getRegion()), state.getSequence());
			if (!t.getReference().getStrand().equals(state.getSequenceRegion().getReference().getStrand()))
				SequenceUtils.getDnaReverseComplementInplace(readSeq);
		}
		else {
			readSeq = genomic.getSequence(t).toString().toUpperCase().toCharArray();
		}

		
		AlignedReadsDataFactory fac = new AlignedReadsDataFactory(t.getData().getNumConditions());
		fac.start();
		

		SubreadsAlignedReadsData data = t.getData();
		for (int d=0; d<data.getDistinctSequences(); d++) {
			
			if (data.hasNonzeroInformation()) {
				int[] nz = data.getNonzeroCountIndicesForDistinct(d);
				for (int i=0; i<nz.length; i++) {
					int count = data.getNonzeroCountValueInt(d, i, ReadCountMode.All);
					for (int c=0; c<count; c++)
						singleResimulate(readSeq,nz[i],t,d,state,fac);
				}
			} else {
				for (int i=0; i<data.getNumConditions(); i++) {
					int count = data.getCount(d, i);
					for (int c=0; c<count; c++)
						singleResimulate(readSeq,i,t,d,state,fac);
				}
			}
			
		}
		if (fac.getDistinctSequences()==0) return null;
		
		fac.makeDistinct();
		ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> re = new ImmutableReferenceGenomicRegion<>(t.getReference(), t.getRegion(),fac.createSubread());
		return re;
	}


	private void singleResimulate(char[] sequence, int cond, ImmutableReferenceGenomicRegion<SubreadsAlignedReadsData> readRegion, int distinct, ResimulateState state, AlignedReadsDataFactory fac) {
		
		SubreadsAlignedReadsData data = readRegion.getData();
		fac.add(data,distinct, v->{
			if (!v.isMismatch()) return v;
			for (int t=0; t<types.length; t++) {
				char genomic = types[t].getGenomic();
				if (state.getCategory().reversesLabel()) {
					genomic = SequenceUtils.getDnaComplement(genomic);
				}
				if (v.getReferenceSequence().charAt(0)==genomic) return null;
			}
			return v;
		},true,false);
		fac.setCount(cond, 1);
		
		// introduced mismatches cannot be within a deletion!
		ArrayGenomicRegion deletion = new ArrayGenomicRegion();
		for (int v=0; v<data.getVariationCount(distinct); v++)
			if (data.isDeletion(distinct, v)) {
				deletion = deletion.union(new ArrayGenomicRegion(data.getDeletionPos(distinct, v),data.getDeletionPos(distinct, v)+data.getDeletion(distinct, v).length()));
			}
		
		int len = readRegion.getRegion().getTotalLength();
		
		for (int t=0; t<types.length; t++) {
			char genomic = types[t].getGenomic();
			char read = types[t].getRead();
			if (state.getCategory().reversesLabel()) {
				genomic = SequenceUtils.getDnaComplement(genomic);
				read = SequenceUtils.getDnaComplement(read);
			}
			
			Double ntr = ntrs[t][cond][state.getCategory().id()].get(state.getTargets().iterator().next());
			if (ntr==null) ntr = 0.0; // e.g., no4sU
			
			boolean isNew = state.getRnd().nextDouble()<ntr;
			double tbunif = isNew && modelType!=ResimulatorModelType.Binom?state.getRnd().nextDouble():-1;
			
			for (int p=0; p<len; p++) {
				if (sequence[p]==genomic && !deletion.contains(p)) {
					
					int subIndex = data.getSubreadIndexForPosition(distinct, p, len);
					int sub = data.getSubreadId(distinct, subIndex);
					ModelStructure model = models[sub][t][cond];
					
					double mmProb = mmMat.getMismatchFrequencyForCondition(cond, sub, genomic, read);
					
					if (model!=null)  {
						switch (modelType) {
						case Binom: mmProb = model.getBinom().getParameter(isNew?"p.conv":"p.err"); break;
						case TbBinom: 
							double l = model.getTBBinom().getParameter("p.err");
							if (isNew) {
								double u = model.getTBBinom().getParameter("p.mconv");
								double shape = model.getTBBinom().getParameter("shape");
								mmProb = Grand3Model.rtbeta(l,u, Math.exp(shape), Math.exp(-shape), tbunif);
							}
							else {
								mmProb = l;
							}
							break;
						default: throw new RuntimeException();
						}
					}
					
					double pr = state.getRnd().nextDouble();
					if (pr<mmProb) {
						fac.addMismatch(p, genomic, read, false);
					} else { // this is not entirely correct, but should do!
						pr-=mmProb;
						for (char mm : SequenceUtils.nucleotides) {
							if (mm!=genomic && mm!=read) {
								mmProb = mmMat.getMismatchFrequencyForCondition(cond, sub, genomic, mm);
								if (pr<mmProb) {
									fac.addMismatch(p, genomic, mm, false);
									break;
								}
								pr-=mmProb;
							}
						}
					}
					
					
				}
			}
		}
		
	}

	






}
