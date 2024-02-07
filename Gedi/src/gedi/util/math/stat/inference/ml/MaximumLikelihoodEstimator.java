package gedi.util.math.stat.inference.ml;

import java.io.File;
import java.io.IOException;

import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;

public interface MaximumLikelihoodEstimator {

	<D> MaximumLikelihoodParametrization fit(MaximumLikelihoodModel<D> model, MaximumLikelihoodParametrization par, D data, String... fixedNames);

	<D> MaximumLikelihoodParametrization computeProfileLikelihoodPointwiseConfidenceIntervals(
			MaximumLikelihoodModel<D> model, 
			MaximumLikelihoodParametrization par,
			D data, String... fixedNames);

	/**
	 * Write prefix in front of each line (e.g. "condition\t")
	 * @param out
	 * @param prefix
	 * @param points
	 * @param nthreads
	 * @return
	 * @throws IOException
	 */
	<D> MaximumLikelihoodParametrization writeProfileLikelihoods(
			MaximumLikelihoodModel<D> model, 
			MaximumLikelihoodParametrization par, 
			D data, 
			LineWriter out,
			String prefix, String... fixedNames) throws IOException;

	default <D> MaximumLikelihoodParametrization writeProfileLikelihoods(
			MaximumLikelihoodModel<D> model, 
			MaximumLikelihoodParametrization par, 
			D data, 
			File f,
			String fixedNames) throws IOException {
		LineWriter out = new LineOrientedFile(f.getPath()).write();
		out.writef("Parameter\t%s\tdeltaLL\n",EI.wrap(par.getParamNames()).concat("\t"));

		writeProfileLikelihoods(model,par,data,out, "",fixedNames);
		out.close();
		return par;
	}
	
}