package gedi.util.genomic.metagene;

import java.util.ArrayList;
import java.util.function.Consumer;

import org.apache.commons.math3.stat.descriptive.UnivariateStatistic;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import gedi.util.ArrayUtils;
import gedi.util.math.stat.descriptive.MeanVarianceOnline;

public class MetageneDefinition {
	private int width;
	ArrayList<MetageneRangeDefinition> regions = new ArrayList<>();
	Consumer<double[]> normalizer = a->{};
	
	MetageneDefinition(int width) {
		this.width = width;
	}
	
	public MetageneRangeDefinition region(String name, MetageneDataProvider provider) {
		return new MetageneRangeDefinition(name, provider);
	}
	
	public MetageneDefinition normalizeSizeFactor(double sf) {
		normalizer = a->ArrayUtils.mult(a, 1/sf);
		return this;
	}
	
	public MetageneDefinition unitize() {
		normalizer = a->ArrayUtils.normalize(a);
		return this;
	}
	
	public MetageneDefinition standardize() {
		normalizer = a->{
			MeanVarianceOnline mv = new MeanVarianceOnline();
			for (double v : a) mv.add(v);
			double mean = mv.getMean();
			double sd = mv.getStandardDeviation();
			for (int i=0; i<a.length; i++)
				a[i] = (a[i]-mean)/sd;
		};
		return this;
	}
	
	
	public Metagene create() {
		return new Metagene(width, normalizer, regions);
	}
	
	
	public class MetageneRangeDefinition {
		String name;
		double width = 1;
		MetageneRegionProvider reg = MetageneRegionProvider.full();
		MetageneDataProvider provider;
		double where = -1;
		int maxLength = -1;
		double fixedFactor = -1;
		UnivariateStatistic shrinkStat = new Mean();
		
		public MetageneRangeDefinition(String name, MetageneDataProvider provider) {
			this.name = name;
			this.provider = provider;
		}
		
		public MetageneDefinition add() {
			regions.add(this);
			return MetageneDefinition.this;
		}
		
		public MetageneRangeDefinition setShrinkStatistics(UnivariateStatistic shrinkStat) {
			this.shrinkStat = shrinkStat;
			return this;
		}
		
		public MetageneRangeDefinition setWidth(double width) {
			this.width = width;
			return this;
		}
		
		public MetageneRangeDefinition setRegionProvider(MetageneRegionProvider reg) {
			this.reg = reg;
			return this;
		}
		
		public MetageneRangeDefinition setTruncateOrPad(double where, int maxLength) {
			this.where = where;
			this.maxLength = maxLength;
			return this;
		}
		
		public MetageneRangeDefinition setFixedScale(double fixedFactor) {
			this.fixedFactor = fixedFactor;
			return this;
		}
	}
}