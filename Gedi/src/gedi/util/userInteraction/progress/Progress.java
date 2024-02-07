package gedi.util.userInteraction.progress;

import java.util.Locale;
import java.util.function.Supplier;


/**
 * First call init, then set options; after that, update progress and in the end, call finish!
 * @author erhard
 *
 */
public interface Progress {

	/**
	 * May return a new instance!
	 * @return
	 */
	Progress init();
	
	Progress setCount(int count);
	Progress setDescription(CharSequence message);
	Progress setDescription(Supplier<CharSequence> message);
	
	CharSequence getDescription();
	
	default Progress setDescriptionf(String format, Object...args) {
		return setDescription(String.format(Locale.US,format, args));
	}
	Progress setProgress(int count);
	Progress incrementProgress();
	
	void updateView(int index, int total);
	void firstView(int total);
	void lastView(int total);
	
	
	boolean isGoalKnown();
	boolean isRunning();
	ProgressManager getManager();
	void setManager(ProgressManager man);
	
	void finish();
	
}
