package gedi.util.io.text.jhp.display;

import java.util.Set;
import java.util.function.Predicate;

import gedi.app.extension.AbstractPrioritiesExtensionPoint;
import gedi.core.workspace.loader.PreloadInfo;
import gedi.util.ReflectionUtils;

public class DisplayTemplateGeneratorExtensionPoint extends AbstractPrioritiesExtensionPoint<PreloadInfo<?,?>,Class<?>,AutoTemplateGenerator<PreloadInfo<?,?>,?>> {

	private static DisplayTemplateGeneratorExtensionPoint instance;

	public static DisplayTemplateGeneratorExtensionPoint getInstance() {
		if (instance==null) 
			instance = new DisplayTemplateGeneratorExtensionPoint();
		return instance;
	}
	
	protected DisplayTemplateGeneratorExtensionPoint() {
		super((Class)AutoTemplateGenerator.class,(pred,p)->find(pred,p.getItemClass()));
	}

	private static Class find(Predicate<Class<?>> pred, Class<?> itemClass) {
		if (pred.test(itemClass))
			return itemClass;
		for (Class c : ReflectionUtils.getImplementedInterfaces(itemClass)) {
			if (pred.test(c))
				return c;
		}
		return null;
	}
	

}

