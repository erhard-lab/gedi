package gedi.slam.javapipeline;


import java.io.File;

import gedi.core.data.reads.ReadCountMode;
import gedi.util.program.GediParameter;
import gedi.util.program.GediParameterSet;
import gedi.util.program.parametertypes.BooleanParameterType;
import gedi.util.program.parametertypes.DoubleParameterType;
import gedi.util.program.parametertypes.EnumParameterType;
import gedi.util.program.parametertypes.FileParameterType;
import gedi.util.program.parametertypes.IntParameterType;
import gedi.util.program.parametertypes.StringParameterType;

public class SlamCheckParameterSet extends GediParameterSet {

	public GediParameter<Integer> nthreads = new GediParameter<Integer>(this,"nthreads", "The number of threads to use for computations", false, new IntParameterType(), Runtime.getRuntime().availableProcessors());
	
	public GediParameter<String> o = new GediParameter<String>(this,"o", "Output table name", false, new StringParameterType());
	
	public GediParameter<File> out = new GediParameter<File>(this,"output", "Output table", false, new FileParameterType());
	
	public GediParameter<Integer> len= new GediParameter<Integer>(this,"len", "Read length", false, new IntParameterType(),50);
	public GediParameter<Double> u = new GediParameter<Double>(this,"u", "U content", false, new DoubleParameterType(),0.3);
	public GediParameter<Double> conv = new GediParameter<Double>(this,"conv", "Specify the conversion rate for estimation", false, new DoubleParameterType(),-1.0);
	public GediParameter<Double> err = new GediParameter<Double>(this,"err", "Specify the error rate for estimation", false, new DoubleParameterType(),-1.0);
	public GediParameter<Double> simconv = new GediParameter<Double>(this,"simconv", "Specify the conversion rate for simulation", false, new DoubleParameterType(),0.023);
	public GediParameter<Double> simerr = new GediParameter<Double>(this,"simerr", "Specify the error rate for simulation", false, new DoubleParameterType(),4E-4);
	public GediParameter<Integer> minEstimateReads = new GediParameter<Integer>(this,"minEstimateReads", "Specify the minimal number of reads to be used for parameter estimation!", false, new IntParameterType(),10000);
	
	public GediParameter<String> p = new GediParameter<String>(this,"p", "Either a file containing proportions, or as a number", false, new StringParameterType());
	public GediParameter<String> r = new GediParameter<String>(this,"r", "Either a file containing read counts, or as a number", false, new StringParameterType());
	public GediParameter<Integer> genes = new GediParameter<Integer>(this,"genes", "Number of genes (or take from -r)", false, new IntParameterType(),-1);

	public GediParameter<Boolean> param = new GediParameter<Boolean>(this,"param", "Estimate binomial mixture parameters only", false, new BooleanParameterType());

	public GediParameter<Boolean> ci = new GediParameter<Boolean>(this,"ci", "Output CI as well", false, new BooleanParameterType());
	public GediParameter<Boolean> beta = new GediParameter<Boolean>(this,"beta", "Output beta parameters as well", false, new BooleanParameterType());
	
}
