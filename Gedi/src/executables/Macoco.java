package executables;

import gedi.macoco.javapipeline.CountEquivalenceClasses;
import gedi.macoco.javapipeline.EquivalenceClassesEffectiveLengths;
import gedi.macoco.javapipeline.EquivalenceClassesMacoco;
import gedi.macoco.javapipeline.EquivalenceClassesMaximumLikelihood;
import gedi.macoco.javapipeline.MacocoParameterSet;
import gedi.util.program.CommandLineHandler;
import gedi.util.program.GediProgram;

public class Macoco {

	
	public static void main(String[] args) {
		
	
		MacocoParameterSet params = new MacocoParameterSet();
		GediProgram pipeline = GediProgram.create("Macoco",
				new CountEquivalenceClasses(params),
//				new EquivalenceClassesMaximumLikelihood(params)
				new EquivalenceClassesMacoco(params),
				new EquivalenceClassesEffectiveLengths(params)
				);
		GediProgram.run(pipeline, new CommandLineHandler("Macoco","Macoco computes maximal coverage consistent estimates of transcript abundances.",args));
		
	}
}
