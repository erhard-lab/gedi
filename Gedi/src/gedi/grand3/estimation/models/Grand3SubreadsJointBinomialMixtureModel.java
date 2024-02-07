package gedi.grand3.estimation.models;

import gedi.grand3.knmatrix.KNMatrix.KNMatrixElement;
import gedi.util.functions.EI;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodModel;
import gedi.util.math.stat.inference.ml.MaximumLikelihoodParametrization;

public class Grand3SubreadsJointBinomialMixtureModel implements MaximumLikelihoodModel<KNMatrixElement[][]>{

	
	private MaximumLikelihoodModel<KNMatrixElement[]>[] submodels;
	
	public Grand3SubreadsJointBinomialMixtureModel(MaximumLikelihoodModel<KNMatrixElement[]>[] submodels) {
		this.submodels = submodels;
	}

	@Override
	public double logLik(KNMatrixElement[][] data, double[] param) {
		double[] buf = new double[3];
		buf[0] = param[0];
		double re = 0;
		for (int i=0; i<submodels.length; i++) {
			System.arraycopy(param, 1+i*2, buf, 1, buf.length-1);
			re+=submodels[i].logLik(data[i], buf);
		}
		return re;
	}

	@Override
	public MaximumLikelihoodParametrization createPar() {
		int subreads = submodels.length;
		String[] paramNames = new String[1+subreads*2];
		double[] params = new double[1+subreads*2];
		String[] paramFormat = new String[1+subreads*2];
		double[] lowerBounds = new double[1+subreads*2];
		double[] upperBounds = new double[1+subreads*2];
		paramNames[0] = "ntr";
		params[0] = 0.1;
		paramFormat[0] = "%.4f";
		lowerBounds[0] = 0;
		upperBounds[0] = 0;
		
		for (int i=0; i<subreads; i++) {
			paramNames[1+i*2] = "p.err."+i;
			paramNames[2+i*2] = "p.conv."+i;
			params[1+i*2] = 4E-4;
			params[2+i*2] = 0.05;
			paramFormat[1+i*2] = "%.4g";
			paramFormat[2+i*2] = "%.4f";
			lowerBounds[1+i*2] = 0;
			lowerBounds[2+i*2] = 0;
			upperBounds[1+i*2] = 0.01;
			upperBounds[2+i*2] = 0.2;
		}
				
		return new MaximumLikelihoodParametrization("binomSubreadsJoint",paramNames,params,paramFormat,lowerBounds,upperBounds);
	}
	
	public MaximumLikelihoodParametrization createPar(MaximumLikelihoodParametrization[] sub) {
		int subreads = submodels.length;
		String[] paramNames = new String[1+subreads*2];
		double[] params = new double[1+subreads*2];
		String[] paramFormat = new String[1+subreads*2];
		double[] lowerBounds = new double[1+subreads*2];
		double[] upperBounds = new double[1+subreads*2];
		paramNames[0] = "ntr";
		params[0] = EI.wrap(sub).mapToDouble(m->m.getParameter("ntr")).mean();
		paramFormat[0] = "%.4f";
		lowerBounds[0] = EI.wrap(sub).mapToDouble(m->m.getLowerBounds("ntr")).max();
		upperBounds[0] = EI.wrap(sub).mapToDouble(m->m.getUpperBounds("ntr")).min();
		
		for (int i=0; i<subreads; i++) {
			paramNames[1+i*2] = "p.err."+i;
			paramNames[2+i*2] = "p.conv."+i;
			params[1+i*2] = sub[i].getParameter(1);
			params[2+i*2] = sub[i].getParameter(2);
			paramFormat[1+i*2] = "%.4g";
			paramFormat[2+i*2] = "%.4f";
			lowerBounds[1+i*2] = sub[i].getLowerBounds(1);
			lowerBounds[2+i*2] = sub[i].getLowerBounds(2);
			upperBounds[1+i*2] = sub[i].getUpperBounds(1);
			upperBounds[2+i*2] = sub[i].getUpperBounds(2);
		}
				
		return new MaximumLikelihoodParametrization("binomSubreadsJoint",paramNames,params,paramFormat,lowerBounds,upperBounds);
	}

	public void distributeTo(MaximumLikelihoodParametrization joint, MaximumLikelihoodParametrization[] para) {
		for (int i=0; i<para.length; i++) {
			para[i].setParameter("ntr", joint.getParameter("ntr"));
			para[i].setParameter("p.err", joint.getParameter("p.err."+i));
			para[i].setParameter("p.conv", joint.getParameter("p.conv."+i));
			if (joint.hasConfidences()) {
				para[i].setLowerConfidence("ntr", joint.getLowerConfidence("ntr"));
				para[i].setLowerConfidence("p.err", joint.getLowerConfidence("p.err."+i));
				para[i].setLowerConfidence("p.conv", joint.getLowerConfidence("p.conv."+i));
				para[i].setUpperConfidence("ntr", joint.getUpperConfidence("ntr"));
				para[i].setUpperConfidence("p.err", joint.getUpperConfidence("p.err."+i));
				para[i].setUpperConfidence("p.conv", joint.getUpperConfidence("p.conv."+i));
			}
		}
	}

}
