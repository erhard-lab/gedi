package gedi.util.math.stat.inference.ml;

import java.util.function.ToDoubleFunction;

public class MaximumLikelihoodModelPriorDecorator<D> implements MaximumLikelihoodModel<D> {

	private ToDoubleFunction<double[]> prior;
	private MaximumLikelihoodModel<D> model;

	public MaximumLikelihoodModelPriorDecorator(MaximumLikelihoodModel<D> model, ToDoubleFunction<double[]> prior) {
		this.prior = prior;
		this.model = model;
	}

	@Override
	public double logLik(D data, double[] param) {
		return prior.applyAsDouble(param)+model.logLik(data, param);
	}

	@Override
	public MaximumLikelihoodParametrization createPar() {
		return model.createPar();
	}
	
	
	
	
}
