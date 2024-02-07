package gedi.core.region.feature.features;

import gedi.core.data.annotation.Transcript;
import gedi.core.region.GenomicRegion;
import gedi.core.region.GenomicRegionPosition;
import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.core.region.ReferenceGenomicRegion;
import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.util.ArrayUtils;
import gedi.util.datastructure.array.NumericArray;

import java.util.Iterator;
import java.util.Set;
import java.util.function.UnaryOperator;


@GenomicRegionFeatureDescription(fromType=ReferenceGenomicRegion.class,toType=Integer.class)
public class RelativePosition extends AbstractFeature<Integer> {

	
	
	private GenomicRegionPosition readPosition = GenomicRegionPosition.FivePrime;
	private int readOffset = 0;
	private double binSize = 0.01;
	
	private boolean forceFullRegion = false;
	private boolean forceCoding = false;
	private double[] medianSizes = null;
	
	public RelativePosition() {
		minInputs = maxInputs = 1;
	}

	public void setReadOffset(int readOffset) {
		this.readOffset = readOffset;
	}
	
	public void setReadPosition(GenomicRegionPosition readPosition) {
		this.readPosition = readPosition;
	}

	public void setBins(int bins) {
		binSize = 1.0/bins;
	}
	
	@Override
	public GenomicRegionFeature<Integer> copy() {
		RelativePosition re = new RelativePosition();
		re.copyProperties(this);
		re.readOffset = readOffset;
		re.readPosition = readPosition;
		re.binSize = binSize;
		re.forceCoding = forceCoding;
		re.forceFullRegion = forceFullRegion;
		return re;
	}
	
	/**
	 * Even for transcripts, force the full region (otherwise the 5'utr, cds and 3'utr are scaled independently)
	 * @param forceFullRegion
	 */
	public void setForceFullRegion(boolean forceFullRegion) {
		this.forceFullRegion = forceFullRegion;
	}
	
	public void setForceCoding(boolean forceCoding) {
		this.forceCoding = forceCoding;
	}
	
	@Override
	public void begin() {
		super.begin();
		if (!forceFullRegion) {
			Iterator<? extends ReferenceGenomicRegion<?>> uit = program.getFeature(getInputName(0)).getUniverse();
			if (uit!=null) {
				medianSizes = new double[3];
				int[] count = new int[3];
				while (uit.hasNext()) {
					ReferenceGenomicRegion<?> n = uit.next();
					if (n.getData() instanceof Transcript) {
						Transcript t = (Transcript) n.getData();
						if (t.isCoding()) {
							addIfNonZero(medianSizes,count,0,t.get5Utr(n.getReference(), n.getRegion()).getTotalLength());
							addIfNonZero(medianSizes,count,1,t.getCds(n.getReference(), n.getRegion()).getTotalLength());
							addIfNonZero(medianSizes,count,2,t.get3Utr(n.getReference(), n.getRegion()).getTotalLength());
						}
					} else 
						break;
				}
				
				if (ArrayUtils.max(count)==0) 
					medianSizes = null;
				else 
					for (int i=0; i<3; i++) 
						medianSizes[i] /= count[i];
			}
		}
		
	}
	

	private void addIfNonZero(double[] a, int[] c, int i,
			int l) {
		if (l>0) {
			a[i]+=l;
			c[i]++;
		}
	}
	
	
	private double currentFac = 1;
	@Override
	public UnaryOperator<NumericArray> getCountAdapter() {
		return na->na.applyInPlace(e->e*currentFac);
	}

	@Override
	protected void accept_internal(Set<Integer> values) {
	
		currentFac = 1;
		
		Set<ReferenceGenomicRegion<?>> inputs = getInput(0);
		
		for (ReferenceGenomicRegion<?> rgr : inputs) {
		
			if (!readPosition.isValidInput(referenceRegion))
				continue;
			
			int r = readPosition.position(referenceRegion,readOffset);
			
			if (!rgr.getRegion().contains(r))
				continue;
			
			if (!forceFullRegion && medianSizes==null) {
				System.err.println("Warning: Input universe could not be iterated! Setting ForceFullRegion!");
				forceFullRegion = true;
			}

			
			if (!forceFullRegion) {
				
				
				Transcript t = ((Transcript) rgr.getData()); 
				if (t.isCoding() && t.get5Utr(rgr.getReference(), rgr.getRegion()).getTotalLength()>0 && t.get3Utr(rgr.getReference(), rgr.getRegion()).getTotalLength()>0) {
					
					for (int p=0; p<3; p++) {
						GenomicRegion part = t.getPart(rgr.getReference(), rgr.getRegion(), p);
						if (part.contains(r)) {
							r = new ImmutableReferenceGenomicRegion(rgr.getReference(), part).induceMaybeOutside(r);
							double rel = (double)r/part.getTotalLength();
							
							// scale rel to median size of this part
							rel = rel/medianSizes[1]*medianSizes[p];
							if (p==0) rel = -medianSizes[0]/medianSizes[1]+rel;
							else if (p==2) rel = 1+rel;
							
							int bin = (int) (rel/binSize);
							values.add(bin);
							
							currentFac = (medianSizes[p]/medianSizes[1])/(part.getTotalLength()/(double)t.getCds(rgr.getReference(), rgr.getRegion()).getTotalLength());
						}
					}
				}
				
			} else {
				
				Transcript t = ((Transcript) rgr.getData()); 
				
				if (!forceCoding || t.isCoding()) {
					r = rgr.induceMaybeOutside(r);
					double rel = (double)r/rgr.getRegion().getTotalLength();
					
					int bin = (int) (rel/binSize);
//					if (bin==81) {
//						sum+=((AlignedReadsData)referenceRegion.getData()).getSumCount(0, true);
//						System.out.println(sum+"\t"+referenceRegion+"\t"+rgr);
//					}
					values.add(bin);
				}
			}
			
		}
		
	}
	
//	public boolean dependsOnData() {
//		return true;
//	}
//
//	double sum = 0;
}

