package executables;

import java.io.File;
import java.io.IOException;

import gedi.slam.GeneData;
import gedi.slam.OptimNumericalIntegrationProportion;
import gedi.slam.ReadData;
import gedi.slam.SlamCollector;
import gedi.slam.SlamEstimationResult;
import gedi.slam.SlamParameterEstimation;
import gedi.slam.javapipeline.SlamCheckParameterSet;
import gedi.util.StringUtils;
import gedi.util.datastructure.array.NumericArray;
import gedi.util.datastructure.array.sparse.AutoSparseDenseDoubleArrayCollector;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.math.stat.NumericSample;
import gedi.util.math.stat.RandomNumbers;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;
import gedi.util.program.GediProgramContext;
import jdistlib.Beta;
import jdistlib.Binomial;
import jdistlib.Normal;
import jdistlib.Uniform;
import jdistlib.rng.MersenneTwister;

public class SlamCheck {

	
	public static void main(String[] args) {
		
		SlamCheckParameterSet params = new SlamCheckParameterSet();
		GediProgram pipeline = GediProgram.create("SlamCheck",
				new SlamChecker(params)
				);
		GediProgram.run(pipeline, new CommandLineHandler("SlamCheck","SlamCheck tests the estimation of p.",args));
		
	}
	
	
	private static class SlamChecker extends GediProgram {
		

		public SlamChecker(SlamCheckParameterSet params) {
			addInput(params.nthreads);
			addInput(params.conv);
			addInput(params.err);
			addInput(params.simconv);
			addInput(params.simerr);
			addInput(params.len);
			addInput(params.u);
			addInput(params.o);
			addInput(params.minEstimateReads);
			addInput(params.param);
			addInput(params.r);
			addInput(params.p);
			addInput(params.genes);
			addInput(params.ci);
			addInput(params.beta);
			
			addOutput(params.out);
			
		}
		
		
		
		public String execute(GediProgramContext context) throws IOException {
			
			int nthreads = getIntParameter(0);
			double conv = getDoubleParameter(1);
			double err = getDoubleParameter(2);
			double sconv = getDoubleParameter(3);
			double serr = getDoubleParameter(4);
			int len = getIntParameter(5);
			double u = getDoubleParameter(6);
			int minEstimateReads = getIntParameter(8);
			boolean paramOnly = getBooleanParameter(9);
			String rr = getParameter(10);
			String pp = getParameter(11);
			int genes = getParameter(12);
			boolean ci = getParameter(13);
			boolean beta = getParameter(14);
			
			String output = getParameter(7);
			
			context.getLog().info("SConv="+sconv);
			context.getLog().info("SErr="+serr);
			context.getLog().info("Len="+len);
			context.getLog().info("U="+u);
			
			int[] readDistr = new File(rr).exists()?EI.lines(rr).mapToInt(Integer::parseInt).toIntArray():new int[] {Integer.parseInt(rr)};
			double[] pDistr = new File(pp).exists()?EI.lines(pp).mapToDouble(Double::parseDouble).toDoubleArray():new double[] {Double.parseDouble(pp)};
			
			if (genes<0) genes = readDistr.length;
			
			int[] reads = EI.seq(0, genes).mapToInt(i->readDistr[RandomNumbers.getGlobal().getUnif(0, readDistr.length)]).toIntArray();
			double[] ps = EI.seq(0, genes).mapToDouble(i->pDistr[RandomNumbers.getGlobal().getUnif(0, pDistr.length)]).toDoubleArray();
			
			
			
			SimulationData[] sim = EI.seq(0, genes).progress(context.getProgress(), genes, i->"Simulate "+i).map(ind->{
//				RandomNumbers rnd = new RandomNumbers();
				MersenneTwister mer = new MersenneTwister();
				
				double p = ps[ind];
				
				ReadData[] rd = new ReadData[reads[ind]];
				for (int r=0; r<rd.length; r++) {
					int total = (int)Binomial.random(len, u, mer);//rnd.getBinom(len, u);
					boolean isnew = Uniform.random(0, 1, mer)<p;//rnd.getUnif()<p;
					int conversions = (int)Binomial.random(total, isnew?sconv:serr, mer);//rnd.getBinom(total, isnew?sconv:serr);
					rd[r] = new ReadData(total, conversions, AutoSparseDenseDoubleArrayCollector.wrap(1.0));
				}
				rd = SlamCollector.collapse(rd);
				
				return new SimulationData(new GeneData("Test", rd),p,reads[ind]);
					
			}).toArray(SimulationData.class);

			if (err<0) 
				err=serr+Normal.random(0, 1E-4, new MersenneTwister());
			
			if (conv<0) {
				double[] param = new SlamParameterEstimation(true).setMinEstimateReads(minEstimateReads)
						.fromGeneData(EI.wrap(sim).map(s->s.gd).toArray(GeneData.class),new double[]{err});
				
				if (conv<0) conv = param[0];
			}
			
			
			context.getLog().info("Conv="+conv);
			context.getLog().info("Err="+err);
			
			boolean ex = new File(output).exists();
			LineWriter writer = new LineOrientedFile(output).append();
			long hash = System.currentTimeMillis();
			
			if (paramOnly) {
				if (!ex)
					writer.write("True conv\tTrue err\tConv\tErr\n");
				writer.writef("%.5g\t%.5g\t%.5g\t%.5g\n", sconv,serr,conv,err);
			}
			else {
				if (!ex) {
					if (beta)
						writer.write("Hash\tTrue conv\tTrue err\tConv\tErr\tReads\ttruth\tinferred\talpha\tbeta\tposterior\n");
					else if (ci)
						writer.write("Hash\tTrue conv\tTrue err\tConv\tErr\tReads\ttruth\tlower\tinferred\tupper\tposterior\n");
					else
						writer.write("Hash\tTrue conv\tTrue err\tConv\tErr\tReads\ttruth\tinferred\tposterior\n");
					
				}
				
				OptimNumericalIntegrationProportion opt = new OptimNumericalIntegrationProportion(1, 1, err, conv,0,0, 0.05, 0.95, true);
				double uconv = conv;
				double uerr = err;
				
				EI.wrap(sim).progress(context.getProgress(), sim.length, i->"Estimate").parallelized(nthreads, 10, ei->ei.map(s->{
					SlamEstimationResult re = opt.infer(s.gd.getConversionVector(0), s.gd.getTotalVector(0), s.gd.getWeightVector(0), new int[0], new int[0], new double[0], context.getLog());
					if (beta)
						return String.format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
								hash,sconv,serr,uconv,uerr,
								s.reads,
								s.p,
								re.getMap(),re.getAlpha(),re.getBeta(),
								Beta.cumulative(s.p, re.getAlpha(), re.getBeta(), true, false));
					else if (ci)
						return String.format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
								hash,sconv,serr,uconv,uerr,
								s.reads,
								s.p,
								re.getLower(),re.getMap(),re.getUpper(),
								Beta.cumulative(s.p, re.getAlpha(), re.getBeta(), true, false));
					else
						return String.format("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%d\t%.4f\t%.4f\t%.4f",
								hash,sconv,serr,uconv,uerr,
								s.reads,
								s.p,
								re.getMap(),
								Beta.cumulative(s.p, re.getAlpha(), re.getBeta(), true, false));
						
				})).print(writer);
			}
			writer.close();
			
			
			return null;
		}


	}
	
	private static class SimulationData {
		GeneData gd;
		double p;
		int reads;
		public SimulationData(GeneData gd, double p, int reads) {
			super();
			this.gd = gd;
			this.p = p;
			this.reads = reads;
		}
		
	}
}
