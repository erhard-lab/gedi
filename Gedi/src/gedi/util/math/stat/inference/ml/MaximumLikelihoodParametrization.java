package gedi.util.math.stat.inference.ml;

import java.io.IOException;

import gedi.util.ArrayUtils;
import gedi.util.GeneralUtils;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;
import gedi.util.io.text.HeaderLine;

public class MaximumLikelihoodParametrization implements BinarySerializable {

	
	String name;
	String[] paramNames;
	double[] params;
	String[] paramFormat;
	double[] lowerBounds;
	double[] upperBounds;

	double logLik;
	String convergence;
	
	double[] lowerConfidence;
	double[] upperConfidence;
	
	
	public MaximumLikelihoodParametrization(String name, String[] paramNames, double[] params, String[] format, double[] lowerBounds, double[] upperBounds) {
		this.name = name;
		this.paramNames = paramNames;
		this.params = params;
		this.paramFormat = format;
		this.lowerBounds = lowerBounds;
		this.upperBounds = upperBounds;
	}
	
	/**
	 * Deep copy only for params and confidence intervals!
	 * @return
	 */
	public MaximumLikelihoodParametrization copy() {
		MaximumLikelihoodParametrization re = new MaximumLikelihoodParametrization(name,paramNames,params.clone(),paramFormat,lowerBounds,upperBounds);
		re.logLik = logLik;
		re.convergence = convergence;
		re.lowerConfidence = lowerConfidence!=null?lowerConfidence.clone():null;
		re.upperConfidence = upperConfidence!=null?upperConfidence.clone():null;
		return re;
	}
	
	/**
	 * True, if all parameters are in the given tolerance (relative to the value of this object)
	 * @param o
	 * @param reltol
	 * @return
	 */
	public boolean isSimilar(MaximumLikelihoodParametrization o, double reltol) {
		for (int i=0; i<params.length; i++)
			if (Math.abs(params[i]-o.params[i])>params[i]*reltol)
				return false;
		return true;
	}
	
	public String getName() {
		return name;
	}
	
	public String[] getParamFormat() {
		return paramFormat;
	}
	
	public String format(int i) {
		return String.format(paramFormat[i], params[i]);
	}
		
	@Override
	public void serialize(BinaryWriter out) throws IOException {
		for (int i=0; i<params.length; i++)
			out.putDouble(params[i]);
		out.putDouble(logLik);
		out.putString(convergence);
		if (lowerConfidence==null) 
			out.putByte(0);
		else {
			out.putByte(1);
			for (int i=0; i<lowerConfidence.length; i++)
				out.putDouble(lowerConfidence[i]);
			for (int i=0; i<upperConfidence.length; i++)
				out.putDouble(upperConfidence[i]);
		}
	}
	
	@Override
	public void deserialize(BinaryReader in) throws IOException {
		for (int i=0; i<params.length; i++)
			params[i] = in.getDouble();
		logLik = in.getDouble();
		convergence = in.getString();
		if (in.getByte()==1) {
			lowerConfidence = new double[params.length];
			upperConfidence = new double[params.length];
			for (int i=0; i<params.length; i++)
				lowerConfidence[i]=in.getDouble();
			for (int i=0; i<params.length; i++)
				upperConfidence[i]=in.getDouble();
		}
	}
	
	public double getLogLik() {
		return logLik;
	}
	
	public int getNumParams() {
		return params.length;
	}
	
	public boolean hasParameter(String name) {
		return ArrayUtils.linearSearch(paramNames, name)!=-1;
	}
		
		
	public int getParameterIndex(String name) {
		int idx = ArrayUtils.linearSearch(paramNames, name);
		if (idx==-1) throw new IndexOutOfBoundsException();
		return idx;
	}
	
	public double getParameter(String name) {
		return params[getParameterIndex(name)];
	}
	
	public double getParameter(int index) {
		return params[index];
	}
	
	public boolean isInBounds() {
		for (int i=0; i<getNumParams(); i++)
			if (getParameter(i)<=getLowerBounds(i) || getParameter(i)>=getUpperBounds(i))
				return false;
		return true;
	}
	
	MaximumLikelihoodParametrization setConfidenceToBounds() {
		for (int i=0; i<getNumParams(); i++) {
			setLowerConfidence(i, getLowerBounds(i));
			setUpperConfidence(i, getUpperBounds(i));
		}
		return this;
	}

	
	public double getUpperBounds(String name) {
		return upperBounds[getParameterIndex(name)];
	}
	public double getUpperBounds(int index) {
		return upperBounds[index];
	}
	
	public double getLowerBounds(String name) {
		return lowerBounds[getParameterIndex(name)];
	}
	public double getLowerBounds(int index) {
		return lowerBounds[index];
	}
	
	public MaximumLikelihoodParametrization setUpperBounds(String name, double u) {
		upperBounds[getParameterIndex(name)] = u;
		return this;
	}
	public MaximumLikelihoodParametrization setUpperBounds(int index, double u) {
		upperBounds[index] = u;
		return this;
	}
	
	public MaximumLikelihoodParametrization setLowerBounds(String name, double l) {
		lowerBounds[getParameterIndex(name)] = l;
		return this;
	}
	public MaximumLikelihoodParametrization setLowerBounds(int index, double l) {
		lowerBounds[index] = l;
		return this;
	}

	public boolean hasConfidences() {
		return lowerConfidence!=null;
	}

	
	public double getUpperConfidence(String name) {
		return upperConfidence[getParameterIndex(name)];
	}
	public double getUpperConfidence(int index) {
		return upperConfidence[index];
	}
	
	public double getLowerConfidence(String name) {
		return lowerConfidence[getParameterIndex(name)];
	}
	public double getLowerConfidence(int index) {
		return lowerConfidence[index];
	}
	
	public MaximumLikelihoodParametrization setUpperConfidence(String name, double l) {
		if (Double.isNaN(l)) l = Double.POSITIVE_INFINITY;
		upperConfidence[getParameterIndex(name)] = l;
		return this;
	}
	public MaximumLikelihoodParametrization setUpperConfidence(int index, double u) {
		if (Double.isNaN(u)) u = Double.POSITIVE_INFINITY;
		upperConfidence[index] = u;
		return this;
	}
	
	public MaximumLikelihoodParametrization setLowerConfidence(String name, double l) {
		if (Double.isNaN(l)) l = Double.NEGATIVE_INFINITY;
		lowerConfidence[getParameterIndex(name)] = l;
		return this;
	}
	public MaximumLikelihoodParametrization setLowerConfidence(int index, double u) {
		if (Double.isNaN(u)) u = Double.NEGATIVE_INFINITY;
		lowerConfidence[index] = u;
		return this;
	}

	
	public MaximumLikelihoodParametrization setParameter(String name, double value) {
		params[getParameterIndex(name)] = value;
		return this;
	}
	
	public MaximumLikelihoodParametrization setParameter(int index, double value) {
		params[index] = value;
		return this;
	}
	
	
	public double[] getParameter() {
		return params;
	}
	
	public String[] getParamNames() {
		return paramNames;
	}
	
	public String getConvergence() {
		return convergence;
	}

	public boolean isConverged() {
		return convergence==null || convergence.length()==0;
	}
	
	public MaximumLikelihoodParametrization addConvergence(String convergence) {
		if (isConverged()) this.convergence = convergence;
		else this.convergence = this.convergence+"; "+convergence;
		return this;
	}

	public void resetConvergence() {
		this.convergence = null;
	}
	

	public String getTabHeader(String prefix) {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append("\t");
			sb.append(prefix).append(paramNames[i]);
		}
		return sb.toString();
	}

	
	public String toTabString() {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append("\t");
			sb.append(String.format(paramFormat[i], params[i]));
		}
		return sb.toString();
	}
	
	public String toTabStringLower() {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append("\t");
			sb.append(String.format(paramFormat[i], lowerConfidence[i]));
		}
		return sb.toString();
	}

	public String toTabStringUpper() {
		StringBuilder sb = new StringBuilder();
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append("\t");
			sb.append(String.format(paramFormat[i], upperConfidence[i]));
		}
		return sb.toString();
	}


	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(name).append("\n");
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append("\n ");
			else sb.append(" ");
			sb.append(paramNames[i])
				.append("=")
				.append(String.format(paramFormat[i], params[i]));
			if (lowerConfidence!=null)
				sb.append(" [").append(String.format(paramFormat[i],lowerConfidence[i]))
					.append(",").append(String.format(paramFormat[i],upperConfidence[i]))
					.append("]");
		}
		sb.append("\nLog likelihood: ").append(String.format("%.2f", logLik));
		if (convergence!=null && convergence.length()>0)
			sb.append("\n").append(convergence);
		return sb.toString();
	}
	
	public String toCreateString(boolean named) {
		StringBuilder sb = new StringBuilder();
		sb.append(getName()).append("(");
		for (int i=0; i<params.length; i++) {
			if (i>0) sb.append(",");
			if (named) sb.append(paramNames[i]).append("=");
			sb.append(String.format(paramFormat[i], params[i]));
		}
		sb.append(")");
		return sb.toString();
	}

	public void parse(String[] a, HeaderLine h, String name) {
		logLik = Double.parseDouble(a[h.get(name+" log likelihood")]);
		convergence="";
		for (int i=0; i<getParamNames().length; i++) 
			params[i] = Double.parseDouble(a[h.get(name+" "+paramNames[i])]);
		
		if (h.hasField("Lower "+name+" "+paramNames[0])) {
			lowerConfidence = new double[params.length];
			upperConfidence = new double[params.length];
			
			for (int i=0; i<getParamNames().length; i++) 
				lowerConfidence[i] = Double.parseDouble(a[h.get("Lower "+name+" "+paramNames[i])]);
			for (int i=0; i<getParamNames().length; i++) 
				upperConfidence[i] = Double.parseDouble(a[h.get("Upper "+name+" "+paramNames[i])]);
		}
		
	}

	public MaximumLikelihoodParametrization integrateOptimization(MaximumLikelihoodParametrization other) {
		this.convergence = GeneralUtils.orDefault(this.convergence, "")+GeneralUtils.orDefault(other.convergence, "");
		this.logLik += other.logLik;
		return this;
	}

	
	
	
	
}
