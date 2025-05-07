package gedi.util.functions;

import java.util.ArrayList;

@SuppressWarnings({"unchecked","rawtypes"})
public class CompoundParallelizedState implements ParallelizedState<CompoundParallelizedState> {

	private ArrayList<ParallelizedState> states = new ArrayList<ParallelizedState>();
	
	public CompoundParallelizedState() {}
	public CompoundParallelizedState(ParallelizedState...states) {
		for (ParallelizedState s : states)
			add(s);
	}
	
	
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
			re.add(s==null?null:s.spawn(index));
		return re;
	}

	@Override
	public void integrate(CompoundParallelizedState other) {
		for (int i=0; i<states.size(); i++)
			if (states.get(i)==null)
				states.set(i, other.states.get(i));
			else if (other.states.get(i)!=null)
				states.get(i).integrate(other.states.get(i));
	}

	public int size() {
		return states.size();
	}

}
