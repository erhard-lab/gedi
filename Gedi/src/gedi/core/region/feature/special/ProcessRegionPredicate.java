package gedi.core.region.feature.special;

import java.util.Set;
import java.util.function.Predicate;

import javax.script.ScriptException;

import gedi.core.region.feature.GenomicRegionFeature;
import gedi.core.region.feature.GenomicRegionFeatureDescription;
import gedi.core.region.feature.features.AbstractFeature;
import gedi.util.nashorn.JSPredicate;

@GenomicRegionFeatureDescription(toType=Boolean.class)
public class ProcessRegionPredicate extends AbstractFeature<Boolean> {

	
	private Predicate<Set[]> pred;
	
	private ProcessRegionPredicate() {
		minValues = maxValues = 1;
	}
	
	public ProcessRegionPredicate(Predicate<Set[]> pred){
		this.pred = pred;
	}
	public ProcessRegionPredicate(String js) throws ScriptException {
		minValues = maxValues = 1;
		StringBuilder code = new StringBuilder();
		code.append("function(i) {\n");
		if (js.contains(";")) {
			code.append(js);
			code.append("}");
		} else {
			code.append("return "+js+";\n}");
		}
		pred = new JSPredicate(code.toString());
	}
	
	@Override
	public boolean dependsOnData() {
		return true;
	}

	@Override
	protected void accept_internal(Set<Boolean> values) {
		values.add(pred.test(inputs));
	}

	@Override
	public GenomicRegionFeature<Boolean> copy() {
		ProcessRegionPredicate re = new ProcessRegionPredicate();
		re.copyProperties(this);
		re.pred = pred;
		return re;
	}
	
}

