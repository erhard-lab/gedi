package gedi.util.io.text.jhp;

import java.util.function.BiConsumer;
import java.util.function.Consumer;

@FunctionalInterface
/**
 * T is the class itself. The array contains either only this, or this and several other template generators
 * @author erhard
 *
 * @param <T>
 */
public interface TemplateGenerator<T extends TemplateGenerator<T>> extends BiConsumer<TemplateEngine,T[]> {
	
}
