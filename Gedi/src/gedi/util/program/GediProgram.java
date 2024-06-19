package gedi.util.program;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.BooleanSupplier;
import java.util.logging.Logger;

import gedi.app.Gedi;
import gedi.util.ArrayUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineOrientedFile;
import gedi.util.io.text.LineWriter;
import gedi.util.job.ExecutionContext;
import gedi.util.job.Job;
import gedi.util.job.PetriNet;
import gedi.util.job.Place;
import gedi.util.job.Transition;
import gedi.util.job.schedule.DefaultPetriNetScheduler;
import gedi.util.mutable.MutableTuple;
import gedi.util.userInteraction.progress.ConsoleProgress;
import gedi.util.userInteraction.progress.NoProgress;
import gedi.util.userInteraction.progress.ProgressManager;

public abstract class GediProgram {
	
	
	protected GediParameterSet parameterSet;
	
	private GediParameterSpec inputSpec = new GediParameterSpec();
	private GediParameterSpec outputSpec = new GediParameterSpec();
	private BooleanSupplier runflag;
	private String name;

	private String changelog = "";
	
	public GediProgram() {
		this.name = getClass().getSimpleName();
	}

	public GediProgram(String name) {
		this.name = name;
	}
	
	public String getChangelog() {
		return changelog;
	}
	
	public void setChangelog(String changelog) {
		this.changelog = changelog;
	}
	
	/**
	 * Runflag: if this is called, the program is only run when the parameter is given. If it is given, it is always! run.
	 * @param parameter
	 */
	protected void setRunFlag(GediParameter<Boolean> parameter) {
		this.runflag = ()->parameter.get().booleanValue();
		addInput(parameter);
	}
	protected void setRunFlag(GediParameter<String> parameter, String runValue) {
		this.runflag = ()->parameter.get().equals(runValue);
		addInput(parameter);
	}
	protected void setRunFlag(BooleanSupplier flag) {
		this.runflag = flag;
	}
	protected <T> void addInput(GediParameter<T> input) {
		inputSpec.add(input);
		if (parameterSet!=null && parameterSet!=input.getParameterSet())
			throw new RuntimeException("Inhomogeneous parameter sets!");
		parameterSet = input.getParameterSet();
	}
	protected <T> void addOutput(GediParameter<T> output) {
		this.outputSpec.add(output);
	}
	
	protected boolean dontRunByUser() {
		return runflag!=null && !runflag.getAsBoolean();
	}
	/**
	 * Does this program want to run even if outputs are already set (default value)?
	 * @return
	 */
	protected boolean wantToRunByUser() {
		return runflag!=null && runflag.getAsBoolean();
	}
	
	protected int getIntParameter(int index) {
		return (Integer)inputSpec.get(index).get();
	}
	
	protected long getLongParameter(int index) {
		return (Long)inputSpec.get(index).get();
	}
	
	protected double getDoubleParameter(int index) {
		return (Double)inputSpec.get(index).get();
	}

	protected boolean getBooleanParameter(int index) {
		return (Boolean)inputSpec.get(index).get();
	}

	protected <T> T getParameter(int index) {
		return (T)inputSpec.get(index).get();
	}
	
	protected <T> ArrayList<T> getParameters(int index) {
		return ( ArrayList<T>)inputSpec.get(index).getList();
	}
	
	protected <T> ArrayList<T> getParameters(String name) {
		return ( ArrayList<T>)inputSpec.get(name).getList();
	}

	protected int getIntParameter(String name) {
		return (Integer)inputSpec.get(name).get();
	}
	
	protected double getDoubleParameter(String name) {
		return (Double)inputSpec.get(name).get();
	}

	protected boolean getBooleanParameter(String name) {
		return (Boolean)inputSpec.get(name).get();
	}

	protected <T> T getParameter(String name) {
		return (T)inputSpec.get(name).get();
	}

	
	protected void setOutput(int index, int value) {
		GediParameter p = outputSpec.get(index);
		p.set(value);
	}
	
	protected void setOutput(int index, double value) {
		GediParameter p = outputSpec.get(index);
		p.set(value);
	}
	
	protected void setOutput(int index, boolean value) {
		GediParameter p = outputSpec.get(index);
		p.set(value);
	}
	
	protected <T> void setOutput(int index, T value) {
		GediParameter p = outputSpec.get(index);
		p.set(value);
	}
	

	public abstract String execute(GediProgramContext context) throws Exception;

	public GediParameterSpec getInputSpec() {
		return inputSpec;
	}

	public GediParameterSpec getOutputSpec() {
		return outputSpec;
	}
	
	public File getOutputFile(int i) {
		return outputSpec.get(i).getFile();
	}

	public LineWriter getOutputWriter(int i) {
		String path = getOutputFile(i).getAbsolutePath();
		return new LineOrientedFile(path).write();
	}

	
	public String getName() {
		return name;
	}
	
	public static void run(GediProgram program, CommandLineHandler cmd) {
		run(program,null,null,cmd);
	}
	public static void run(GediProgram program, GediParameter<File> parameterFile, CommandLineHandler cmd) {
		run(program,parameterFile,null,cmd);
	}
	public static void run(GediProgram program, GediParameter<File> parameterFile, GediParameter<File> runtimeFile, CommandLineHandler cmd) {
		
		try {
			
			Gedi.startup(true);

			String error = cmd.parse(program.getInputSpec(),program.parameterSet);

			if (program.getBooleanParameter(CommandLineHandler.hhh) || program.getBooleanParameter(CommandLineHandler.hh) || program.getBooleanParameter(CommandLineHandler.h) || error!=null) {
				int verbosity = 1;
				if (program.getBooleanParameter(CommandLineHandler.hh)) verbosity=2;
				if (program.getBooleanParameter(CommandLineHandler.hhh)) verbosity=3;
						
				cmd.usage(error,program.getInputSpec(),program.getOutputSpec(),verbosity);
				System.exit(error!=null?1:0);
			}
			
			if (program.getBooleanParameter(CommandLineHandler.changelog)) {
				System.err.println("\n"+program.getChangelog()+"\n");
				System.exit(0);
			}
				
				
			if (program.dontRunByUser())
				return;
			
			if (parameterFile!=null) {
				parameterFile.getFile().getAbsoluteFile().getParentFile().mkdirs();
				program.inputSpec.writeParameterFile(parameterFile.getFile());
			}
			
			program.initParameter(program.parameterSet,program.getInputSpec());
			
			ProgressManager man = new ProgressManager();
			program.getParameter("progress");
			
			GediProgramContext context = new GediProgramContext(
					Gedi.getLog(),
					program.getBooleanParameter(CommandLineHandler.progress)?()->new ConsoleProgress(System.err,man):()->new NoProgress(),
					runtimeFile!=null?runtimeFile.getFile():null,
					program.getBooleanParameter(CommandLineHandler.dry),
					program.parameterSet);
			
			context.writeIntoRunFile("Started %s",program.getName());
			error = program.execute(context);
			context.writeIntoRunFile("Finished %s",program.getName());
			
			if (error!=null) throw new RuntimeException("Could not run program: "+error);
			
		} catch (Throwable e) {
			
			String emsg = StringUtils.createExceptionMessage(e);
			
			cmd.usage(emsg,program.getInputSpec(),program.getOutputSpec(),1);
			
			if (program.getBooleanParameter(CommandLineHandler.D) || ArrayUtils.linearSearch(cmd.getArgs(),"-D")>=0)
				e.printStackTrace();
			System.exit(2);
		}
		
		if (!program.getBooleanParameter(CommandLineHandler.keep) && !program.getBooleanParameter(CommandLineHandler.dry)) {
			for (GediParameter<?> out : program.outputSpec.list)
				if (out.isRemoveFile() && out.getFile().exists()) {
					Logger.getLogger("GEDI").info("Removing temp file "+out.getFile());
					out.getFile().delete();
				}
		}
	}
	
	
	protected void initParameter(GediParameterSet params, GediParameterSpec inputSpec) {
		
	}

	/**
	 * Create a program composed of subprograms.
	 * @return
	 */
	@SuppressWarnings("rawtypes")
	public static GediProgram create(String name, GediProgram...subs1) {
		
		GediProgram[] subs = EI.wrap(subs1).removeNulls().toArray(GediProgram.class);

		HashSet<GediParameter> inputs = new HashSet<>();
		HashSet<GediParameter> outputs = new HashSet<>();
		GediParameterSpec inputSpec = new GediParameterSpec();
		GediParameterSpec outputSpec = new GediParameterSpec();
		
		for (GediProgram s : subs) {
			inputs.addAll(s.inputSpec.list);
			outputs.addAll(s.outputSpec.list);
		}
		HashSet<GediParameter> optionals = EI.wrap(inputs).filter(p->p.isOptional()).set();
		
		EI.wrap(inputs).filter(p->outputs.contains(p)).forEachRemaining(p->p.optional=true);
		EI.wrap(outputs).filter(p->inputs.contains(p)).forEachRemaining(p->p.optional=true);
		
		for (GediProgram s : subs) {
			for (GediParameter in : s.inputSpec.list)
				if (inputs.contains(in))
					inputSpec.addOrMoveToGeneral(s.getName(), in);
			for (GediParameter out : s.outputSpec.list)
				if (outputs.contains(out) & !outputSpec.contains(out))
					outputSpec.add(s.getName(), out);
		}
		
		
		ExecutorService pool = Executors.newCachedThreadPool(new ProgramThreadFactory());
		
		GediProgram re = new GediProgram(name) {
			
			@Override
			protected void initParameter(GediParameterSet params, GediParameterSpec inputSpec) {
				for (GediProgram s : subs) 
					s.initParameter(params, inputSpec);
			}
			
			@Override
			public String execute(GediProgramContext context) throws Exception {
				
				PetriNet pn = new PetriNet();
				HashMap<GediParameter,Place> places = new HashMap<>();
				
				// create a transition for each sub program; for each parameter, connect the corresponding input place
				// if it has a single output, connect the corresponding place; if there are several, connect the transition
				// to a dummy place p and p via dummy transitions to output places
				for (GediProgram s : subs) {
					if (s.dontRunByUser())
						continue;
					GediProgramJob j = new GediProgramJob(s).setDry(context.isDryRun());
					Transition t = pn.createTransition(j);
					for (int i=0; i<s.inputSpec.list.size(); i++) {
						Place p = places.computeIfAbsent(s.inputSpec.list.get(i), x->pn.createPlace(Boolean.class));
						pn.connect(p, t, i);
					}
					if (s.outputSpec.list.size()==1) {
						Place p = places.computeIfAbsent(s.outputSpec.list.get(0), x->pn.createPlace(Boolean.class));
						if (s.outputSpec.list.get(0).get()==null || s.wantToRunByUser()) {
							pn.connect(t,p);
						}else
							j.disable();
					} else if (s.outputSpec.list.size()>1) {
						int dc = 0;
						Place dp = pn.createPlace(Boolean.class);
						pn.connect(t, dp);
						
						for (int i=0; i<s.outputSpec.list.size(); i++) {
							if (s.outputSpec.list.get(i).get()==null || s.wantToRunByUser()) {
								Place p = places.computeIfAbsent(s.outputSpec.list.get(i), x->pn.createPlace(Boolean.class));
								Transition dt = pn.createTransition(new DummyJob());
								pn.connect(dp, dt,0);
								pn.connect(dt, p);
							} else dc++;
						}
						if (dc==s.outputSpec.list.size())
							j.disable();
					}
				}
				
				// create a transition for each place that is not set or optional
				for (GediParameter<?> inp : new ArrayList<>(places.keySet())) {
					if (inp.get()!=null || optionals.contains(inp)) {
						if (places.get(inp)!=null && places.get(inp).getProducer()==null) {
							Transition t = pn.createTransition(new InitJob(inp));
							pn.connect(t, places.computeIfAbsent(inp, x->pn.createPlace(Boolean.class)));
						}
					}
				}
				
				// connect all output places to a finishing transition and place to check whether all outputs were generated
				int ins = (int) EI.wrap(outputs).filter(op->!op.isOptional() && op.get()==null && places.get(op)!=null).count();
				if (ins==0) return "Nothing to do! If you want to re-run clear the output files!"; 
				Transition tfinish = pn.createTransition(new DummyJob(ins)); 
				Place pfinish = pn.createPlace(Boolean.class); 
				pn.connect(tfinish, pfinish);
				int ind = 0;
				for (GediParameter<?> op : outputs)
					if (!op.isOptional() && op.get()==null && places.get(op)!=null)
						pn.connect(places.get(op), tfinish,ind++);

				
				pn.createMissingPlaces();
				
				StringBuilder sb = new StringBuilder();
				HashSet<Place> sources = new HashSet<Place>(pn.getSources());
				for (GediParameter<?> param : places.keySet()) {
					if (sources.contains(places.get(param)))
						sb.append(" ").append(param.getName());
				}
				if (sb.length()>0) return "Mandatory parameters not set:"+sb.toString();
				
				
				pn.prepare();
				
				ExecutionContext econtext = pn.createExecutionContext()
						.newContext(ExecutionContext.UID, String.class)
						.newContext("context", ExecutionContext.class);
				econtext.setContext("context", context);
				DefaultPetriNetScheduler scheduler = new DefaultPetriNetScheduler(econtext, pool);
				scheduler.setRethrowExceptions(true);
//				SingleThreadPetriNetScheduler scheduler = new SingleThreadPetriNetScheduler(econtext);
						
				for (Transition t : pn.getTransitions()) {
					if (t.getJob().isDisabled(econtext))
						econtext.setDisabled(t, true);
				}
				econtext.disableUnneccessary();
				
				HashMap<Place,GediParameter> iplaces = new HashMap<>();
				for (GediParameter<?> param : places.keySet()) 
					iplaces.put(places.get(param), param);
					
				scheduler.run();
				
				return econtext.getToken(pfinish)!=null?null:"Nothing to do!";
			}
		};
		re.inputSpec = inputSpec;
		re.outputSpec = outputSpec;
		re.parameterSet = EI.wrap(subs).map(p->p.parameterSet).unique(true).getUniqueResult("Sub programs must have the same parameterset!", "No subs given!");
		return re;
	}
	


	private static AtomicInteger index = new AtomicInteger(0);
	private static class InitJob implements Job<Boolean> {

		private GediParameter<?> input;
		
		protected int ind;
		
		public InitJob(GediParameter<?> input) {
			this.input = input;
			this.ind = index.getAndIncrement();
		}

		@Override
		public Class[] getInputClasses() {
			return new Class[]{Void.class};
		}

		@Override
		public Class<Boolean> getOutputClass() {
			return Boolean.class;
		}

		@Override
		public boolean isDisabled(ExecutionContext context) {
			return input.get()==null && !input.isOptional();
		}

		@Override
		public Boolean execute(ExecutionContext context, MutableTuple input) {
			return true;
		}

		@Override
		public String getId() {
			return "init"+ind;
		}
	}
	
	private static class ProgramThreadFactory implements ThreadFactory {
        private static final AtomicInteger poolNumber = new AtomicInteger(1);
        private final ThreadGroup group;
        private final AtomicInteger threadNumber = new AtomicInteger(1);
        private final String namePrefix;

        ProgramThreadFactory() {
            SecurityManager s = System.getSecurityManager();
            group = (s != null) ? s.getThreadGroup() :
                                  Thread.currentThread().getThreadGroup();
            namePrefix = "GediProgram-" +
                          poolNumber.getAndIncrement() +
                         "-thread-";
        }

        public Thread newThread(Runnable r) {
            Thread t = new Thread(group, r,
                                  namePrefix + threadNumber.getAndIncrement(),
                                  0);
            t.setDaemon(true);
            t.setPriority(3);
            return t;
        }
    }
}
