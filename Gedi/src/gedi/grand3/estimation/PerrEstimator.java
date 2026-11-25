package gedi.grand3.estimation;

import static gedi.grand3.estimation.MismatchMatrix.genomic;
import static gedi.grand3.estimation.MismatchMatrix.mmIndex;
import static gedi.grand3.estimation.MismatchMatrix.read;

import java.io.File;
import java.io.IOException;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import gedi.grand3.estimation.estimators.ErrorEstimationMethod;
import gedi.grand3.experiment.ExperimentalDesign;
import gedi.grand3.experiment.MetabolicLabel.MetabolicLabelType;
import gedi.util.ArrayUtils;


public class PerrEstimator extends MismatchMatrix {
	
	private static final Median median = new Median();
	private static final Mean mean = new Mean();
	private static final StandardDeviation sd = new StandardDeviation();

	private boolean excludeAntisense = true;
	private double[] manual;
	
	public PerrEstimator(File mmFile, ExperimentalDesign design, String[] subreads, double[] manual) throws IOException {
		super(mmFile,design,subreads);
		this.manual = manual;
	}
	
	public void setExcludeAntisense(boolean excludeAntisense) {
		this.excludeAntisense = excludeAntisense;
	}
	
	
	public double[] estimate(int sample, int subread, MetabolicLabelType type, ErrorEstimationMethod errMethod) {
		int[] no = design.getSamplesNotHaving(type);
		
		if (errMethod.equals(ErrorEstimationMethod.RegressionNo4sU)) {
			if (no.length>0){
				double[] re = regress(mmtab[sample][subread],mmtab[no[0]][subread],type);
				for (int i=1; i<no.length; i++) {
					double[] re2 = regress(mmtab[sample][subread],mmtab[no[i]][subread],type);
					if (re2[1]-re2[0]<re[1]-re[0])
						re = re2;
				}
				if (re[1]-re[0]>0)
					return re;
			}
			throw new RuntimeException("Could not estimate p_err using no4sU samples!");
		}
		
		if (errMethod.equals(ErrorEstimationMethod.SimpleNo4sU)) {
			if (no.length>0){
				double[] rates = new double[no.length];
				for (int i=0; i<rates.length; i++)
					rates[i] = mmtab[no[i]][subread][mmIndex(type.getGenomic(),type.getRead())];
				double mu = mean.evaluate(rates);
				double sem = sd.evaluate(rates)/Math.sqrt(rates.length);
				if(sem<mu*0.1) sem = mu*0.1;
				return new double[] {mu-sem,mu+sem};
			}
			throw new RuntimeException("Could not estimate p_err using no4sU samples!");
		}
		
		if (errMethod.equals(ErrorEstimationMethod.MedianErrorRates)) {
		
			double[] f = flatten(mmtab[sample][subread],type);
			double med = median.evaluate(f);
			double d = Math.min(ArrayUtils.max(f)-med,med-ArrayUtils.min(f));
			double[] re = new double[]{med-d,med+d};
			
			if (re[1]-re[0]>0)
				return re;

			throw new RuntimeException("Could not estimate p_err using median approach!");
		}
		
		if (errMethod.equals(ErrorEstimationMethod.Manual)) {
			return manual;
		}
		
		throw new RuntimeException("Could not estimate p_err, unknown method!");
	}
	
	private double[] regress(double[] s, double[] on, MetabolicLabelType type) {
		s = flatten(s,type);
		double[] n = flatten(on,type);
		
		SimpleRegression reg = new SimpleRegression(false);
		for (int i=0; i<s.length; i++)
			reg.addData(n[i], s[i]);
		
		double ci = reg.getSlopeConfidenceInterval(0.01);
		double sl = reg.getSlope();
		return new double[] {on[mmIndex(type.getGenomic(),type.getRead())]*(sl-ci),on[mmIndex(type.getGenomic(),type.getRead())]*(sl+ci)};
	}
	
	private double[] flatten(double[] mms, MetabolicLabelType type) {
		double[] re = new double[10];
		int index = 0;
		for (int i=0; i<16; i++) {
			if (genomic(i)!=read(i) && (genomic(i)!=type.getGenomic() || read(i)!=type.getRead()) && (!excludeAntisense || genomic(i)!=type.getGenomicReverse() || read(i)!=type.getReadReverse()))
				re[index++]=mms[i];
		}
		if (index!=10)
			throw new RuntimeException("Fatal exception in p.err estimator ("+type+")");
		return re;
	}
	
}
