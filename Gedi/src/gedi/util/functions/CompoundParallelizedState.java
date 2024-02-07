package gedi.util.functions;

import java.util.ArrayList;

@SuppressWarnings({"unchecked","rawtypes"})
public class CompoundParallelizedState implements ParallelizedState<CompoundParallelizedState> {

	private ArrayList<ParallelizedState> states = new ArrayList<ParallelizedState>();
	
	public <T extends ParallelizedState<T>> CompoundParallelizedState add(T state) {
		states.add(state);
		return this;
	}
	
	public <T extends ParallelizedState<T>> T get(int index) {
		return (T) states.get(index);
	}
	
	@Override
	public CompoundParallelizedState spawn(int index) {
		CompoundParallelizedState re = new CompoundParallelizedState();
		for (ParallelizedState s : states)
			re.add(s.spawn(index));
		return re;
	}

	@Override
	public void integrate(CompoundParallelizedState other) {
		for (int i=0; i<states.size(); i++)
			states.get(i).integrate(other.states.get(i));
	}

	public int size() {
		return states.size();
	}

}
