package gedi.grand3.targets;

import java.util.function.UnaryOperator;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.ExtendedIterator;

/**
 * Decorator
 */
public class TargetCollectionMappedName implements TargetCollection {

	private TargetCollection parent;
	private UnaryOperator<String> mapper;
	
	public TargetCollectionMappedName(TargetCollection parent, UnaryOperator<String> mapper) {
		this.parent = parent;
		this.mapper = mapper;
	}

	@Override
	public ImmutableReferenceGenomicRegion<String> getRegion(String name) {
		throw new RuntimeException("Cannot map name for targetted inference!");
	}

	@Override
	public ExtendedIterator<ImmutableReferenceGenomicRegion<String>> iterateRegions() {
		return parent.iterateRegions().map(r->new ImmutableReferenceGenomicRegion<String>(r.getReference(),r.getRegion(),mapper.apply(r.getData())));
	}

	@Override
	public int getNumRegions() {
		return parent.getNumRegions();
	}

	@Override
	public void classify(ImmutableReferenceGenomicRegion<String> target,
			ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, Strandness strandness,
			boolean isStrandCorrected, Grand3ReadClassified classified) {
		parent.classify(target, read, strandness, isStrandCorrected, classified);
	}

	@Override
	public int getNumCategories() {
		return parent.getNumCategories();
	}

	@Override
	public CompatibilityCategory getCategory(int index) {
		return parent.getCategory(index);
	}

	@Override
	public TargetCollection create(ReadCountMode mode, ReadCountMode overlap) {
		return new TargetCollectionMappedName(parent.create(mode, overlap), mapper);
	}
	
	
	
	
}
