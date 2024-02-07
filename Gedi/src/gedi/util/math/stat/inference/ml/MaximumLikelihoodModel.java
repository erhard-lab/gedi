package gedi.util.math.stat.inference.ml;

public interface MaximumLikelihoodModel<D> {

	double logLik(D data, double[] param);
	MaximumLikelihoodParametrization createPar();
	
	default double logLik(D data, MaximumLikelihoodParametrization param) {
		double re = logLik(data, param.getParameter());
		param.logLik = re;
		return re;
	}
	
}
