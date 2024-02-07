package gedi.grand3.targets;

import java.util.function.Predicate;

import gedi.core.data.reads.AlignedReadsData;
import gedi.core.data.reads.ReadCountMode;
import gedi.core.reference.Strandness;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.functions.EI;
import gedi.util.functions.ExtendedIterator;

/**
 * The mechanics are like this:
 * 1. iterate via {@link #iterateRegions()}; the string contained denotes the target that is supposed to found in this region
 * 2. for each read found in this region: determine compatible targets (how is specified by parameters, e.g. strictly contained etc.)
 * 3. map the target to a name
 * 4. Count it, if the category is appropriate (e.g. exonic only), and the name of a target matches the region name
 * 5. How much to count: depends on parameters mode and overlap
 * 
 * For standard mode (genes and transcripts) this is:
 * Iterate genes, compatible with any transcript, map to gene name
 * 
 * For 10x (identify read clusters first):
 * Iterate clusters, only a single compatible target (the one with biggest overlap, if there are many), assign cluster category
 * @param regions
 * @param targets
 * @param targetToName
 */
public interface TargetCollection {
	
	
	ExtendedIterator<ImmutableReferenceGenomicRegion<String>> iterateRegions();
	int getNumRegions();
	void classify(ImmutableReferenceGenomicRegion<String> target, ImmutableReferenceGenomicRegion<? extends AlignedReadsData> read, Strandness strandness, boolean isStrandCorrected, Grand3ReadClassified classified);
	int getNumCategories();
	CompatibilityCategory getCategory(int index);	
	TargetCollection create(ReadCountMode mode, ReadCountMode overlap);
	
	
	default CompatibilityCategory[] getCategories(Predicate<CompatibilityCategory> predicate) {
		return EI.seq(0, getNumCategories()).map(i->getCategory(i)).filter(predicate).toArray(CompatibilityCategory.class);
	}
	
	/**
	 * The {@link CompatibilityCategory#id()} values exactly reflect the indices herein!
	 * @return
	 */
	default CompatibilityCategory[] getCategories() {
		return EI.seq(0, getNumCategories()).map(i->getCategory(i)).toArray(CompatibilityCategory.class);
	}
	
	
	default void checkValid() {
		for (int i=0; i<getNumCategories(); i++)
			if (getCategories()[i].id()!=i)
				throw new RuntimeException("Fatal error: Category indices do not match!");
	}
	
}
