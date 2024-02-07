package gedi.core.region.feature;

import gedi.core.data.annotation.Transcript;
import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.AlignedReadsDataFactory;
import gedi.core.data.reads.ConditionMappedAlignedReadsData;
import gedi.core.data.reads.ContrastMapping;
import gedi.core.data.reads.SelectDistinctSequenceAlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Chromosome;
import gedi.core.region.ArrayGenomicRegion;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionStorage;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.MutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.feature.index.WriteCoverageRmq;
import gedi.core.region.feature.index.WriteJunctionCit;
import gedi.core.region.feature.output.FeatureListOutput;
import gedi.core.region.feature.output.FeatureStatisticOutput;
import gedi.util.ArrayUtils;
import gedi.util.FileUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.NumericArray.NumericArrayType;
import gedi.util.datastructure.array.functions.NumericArrayFunction;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.Progress;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.function.Consumer;
import java.util.function.Supplier;

public class AlignedReadsDataToFeatureProgram implements Consumer<ReferenceGenomicRegion<AlignedReadsData>> {

	private GenomicRegionFeatureProgram<AlignedReadsData> program;
	private MutableReferenceGenomicRegion<AlignedReadsData> mut = new MutableReferenceGenomicRegion<AlignedReadsData>();
	
	private Progress progress = new NoProgress();

	private ReadCountMode readCountMode = ReadCountMode.Weight;
	private Supplier<CharSequence> descr;
	
	private boolean normalize = false;
	private String normalizeName = "";
	private ContrastMapping contrasts;
	private NumericArray rezSizeFactors; // must be for the contrast mapped conditions!
	private String[] conditions;
	
	public AlignedReadsDataToFeatureProgram(
			GenomicRegionFeatureProgram<AlignedReadsData> program) {
		this.program = program;
		
	}
	
	
	public static AlignedReadsDataToFeatureProgram getSimpleProgram(String statistic, String list, GenomicRegionFeature...features) {
		GenomicRegionFeatureProgram<AlignedReadsData> re = new GenomicRegionFeatureProgram<>();
		
		for (GenomicRegionFeature f : features) {
			if (f.getId()==null)
				f.setId(f.getClass().getSimpleName().replaceAll("Feature$", ""));
			re.add(f);
		}
		
		
		if (statistic!=null)
			re.add(new FeatureStatisticOutput(statistic), EI.wrap(features).map(f->f.getId()).toArray(String.class));
		if (list!=null)
			re.add(new FeatureListOutput(list), EI.wrap(features).map(f->f.getId()).toArray(String.class));
		re.setThreads(2);
		
		return new AlignedReadsDataToFeatureProgram(re);
	}


	public static AlignedReadsDataToFeatureProgram getCoverageProgram(String coverageRmq, String junctionCoverage, boolean normalize, String...conditions) {
		GenomicRegionFeatureProgram<AlignedReadsData> re = new GenomicRegionFeatureProgram<>();
		if (coverageRmq!=null) {
			WriteCoverageRmq f = new WriteCoverageRmq(coverageRmq);
			f.setId("coverage");
			re.add(f);
		}
		if (junctionCoverage!=null) {
			WriteJunctionCit f = new WriteJunctionCit(junctionCoverage);
			f.setId("junction");
			re.add(f);
		}
		re.setThreads(0);
		AlignedReadsDataToFeatureProgram re2 = new AlignedReadsDataToFeatureProgram(re);
		re2.normalize = normalize;
		if (conditions.length>0) 
			re2.conditions = conditions;
		
		
		return re2;
	}
	
	public static AlignedReadsDataToFeatureProgram getCoverageProgramsPerCondition(String coverageRmq, String junctionCoverage, boolean normalize, String...conditions) {
		GenomicRegionFeatureProgram<AlignedReadsData> re = new GenomicRegionFeatureProgram<>();
		if (coverageRmq!=null) {
			for (String c : conditions) {
				WriteCoverageRmq f = new WriteCoverageRmq(FileUtils.getFullNameWithoutExtension(coverageRmq)+"."+c+"."+FileUtils.getExtension(coverageRmq),c);
				f.setId("coverage."+c);
				re.add(f);
			}
		}
		if (junctionCoverage!=null) {
			for (String c : conditions) {
				WriteJunctionCit f = new WriteJunctionCit(FileUtils.getFullNameWithoutExtension(junctionCoverage)+"."+c+"."+FileUtils.getExtension(junctionCoverage),c);
				f.setId("junction."+c);
				re.add(f);
			}
		}
		re.setThreads(0);
		AlignedReadsDataToFeatureProgram re2 = new AlignedReadsDataToFeatureProgram(re);
		re2.normalize = normalize;
		if (conditions.length>0) 
			re2.conditions = conditions;
		
		
		return re2;
	}
	
	public AlignedReadsDataToFeatureProgram setNormalizeName(String normalizeName) {
		this.normalizeName = normalizeName;
		return this;
	}


	public GenomicRegionFeatureProgram<AlignedReadsData> getProgram() {
		return program;
	}
	
	public AlignedReadsDataToFeatureProgram setProgress(Progress progress) {
		this.progress = progress;
		return this;
	}

	public AlignedReadsDataToFeatureProgram setProgress(Progress progress, Supplier<CharSequence> descr) {
		this.progress = progress;
		this.descr = descr;
		return this;
	}

	public AlignedReadsDataToFeatureProgram setConsoleProgress() {
		this.progress = new ConsoleProgress(System.err);
		return this;
	}
	
	public void setReadCountMode(ReadCountMode readCountMode) {
		this.readCountMode = readCountMode;
	}
	
	public ReadCountMode getReadCountMode() {
		return readCountMode;
	}

	public void setSize(int size) {
		this.progress.setCount(size);
	}

	public void test(GenomicRegionStorage<Transcript> anno) {
		Chromosome ref = Chromosome.obtain("chr1-");
		MutableReferenceGenomicRegion<Transcript> rgr = EI.wrap(anno.iterateMutableReferenceGenomicRegions(ref))
				.filter(t->t.getData().isCoding() && t.getData().get3Utr(t.getReference(), t.getRegion()).getTotalLength()>500)
				.skip(100)
				.next();
		GenomicRegion reg = new ImmutableReferenceGenomicRegion(ref, rgr.getData().get3Utr(rgr.getReference(), rgr.getRegion())).map(new ArrayGenomicRegion(10,38));
		MutableReferenceGenomicRegion<AlignedReadsData> read = new MutableReferenceGenomicRegion<AlignedReadsData>();
		read.set(ref, reg, new AlignedReadsDataFactory(1).start().newDistinctSequence().setMultiplicity(1).setCount(0, 1).create());
		
		program.setDataToCounts((ard,b)->ard.getCountsForDistinct(b, 0, readCountMode));//ard.addCount(0, b, true));
		program.begin();
		accept(read);
		program.end();
	}
	
	/**
	 * r itself can be mutable, but the data should not be changed (by reflection or orm), as processing may be concurrent!
	 */
	@Override
	public void accept(ReferenceGenomicRegion<AlignedReadsData> r) {
		
		for (int d=0; d<r.getData().getDistinctSequences(); d++) {
			program.accept(mut.set(r.getReference(),r.getRegion(),new SelectDistinctSequenceAlignedReadsData(r.getData(),d)));
		}
		if (descr!=null)
			progress.setDescription(descr);
		progress.incrementProgress();
	}
	
	public void processStorage(GenomicRegionStorage<AlignedReadsData> storage) {
		program.setLabelsFromStorage(storage);
		program.setDataToCounts((ard,b)->{
			if (contrasts!=null) ard = new ConditionMappedAlignedReadsData(ard, contrasts);
			NumericArray re = ard.getCountsForDistinct(b, 0, readCountMode);
			if (rezSizeFactors!=null) re.mult(rezSizeFactors);
			return re;
		});//ard.addCount(0, b, true));
		
		if (normalize) {
			rezSizeFactors = NumericArray.wrap(storage.getMetaDataTotals(normalizeName));
			double med = rezSizeFactors.evaluate(NumericArrayFunction.Median);
			for (int i=0; i<rezSizeFactors.length(); i++)
				rezSizeFactors.setDouble(i, med/rezSizeFactors.getDouble(i));
		}
		HashSet<String> conditionSet = new HashSet<>(Arrays.asList(storage.getMetaDataConditions()));
		
		if (conditions!=null && contrasts==null) {
			conditionSet.clear();
			HashMap<String, Integer> index = ArrayUtils.createIndexMap(storage.getMetaDataConditions());
			
			NumericArray mappedRezSizeFactors = NumericArray.createMemory(conditions.length, NumericArrayType.Double);
			contrasts = new ContrastMapping();
			for (int i=0; i<conditions.length; i++) {
				conditionSet.add(conditions[i]);
				contrasts.addMapping(index.get(conditions[i]), i, conditions[i]);
				if (rezSizeFactors!=null)
					mappedRezSizeFactors.setDouble(i, rezSizeFactors.getDouble(index.get(conditions[i])));
			}
			if (rezSizeFactors!=null)
				rezSizeFactors = mappedRezSizeFactors;
		}
		
		
		for (WriteCoverageRmq cov : program.getFeatures().castFiltered(WriteCoverageRmq.class).loop()) {
			if (cov!=null)
				try {
					FileUtils.writeAllText(
							DynamicObject.from("conditions", DynamicObject.from(
										EI.wrap(storage.getMetaData().getEntry("conditions").asArray())
											.filter(dyn->(cov.getCondition()==null || cov.getCondition().equals(dyn.getEntry("name").asString())) && conditionSet.contains(dyn.getEntry("name").asString()))
											.iff(normalize, ei->ei.map(dyn->dyn.removeEntries(en->en.startsWith("total"))))
											.toArray()
										)).toJson(),
							new File(cov.getPath()+".metadata.json")
							);
				} catch (IOException e) {
					throw new RuntimeException("Could not write metadata!");
				}
		}

		
		program.begin();
		progress.init();
		progress.setCount((int) storage.size());
		storage.ei().forEachRemaining(this);
		progress.finish();
		program.end();
		

	}
	
	public void process(Iterator<? extends ReferenceGenomicRegion<AlignedReadsData>> it) {
		process(it,-1);
	}
	
	public void process(Iterator<? extends ReferenceGenomicRegion<AlignedReadsData>> it, int count) {
		program.setDataToCounts((ard,b)->ard.getCountsForDistinct(b, 0, readCountMode));//ard.addCount(0, b, true));
		
		program.begin();
		progress.init();
		if (count>0)
			progress.setCount(count);
		while (it.hasNext())
			accept(it.next());
		progress.finish();
		program.end();

	}

}
