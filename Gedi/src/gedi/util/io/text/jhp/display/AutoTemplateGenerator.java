package gedi.util.io.text.jhp.display;

import java.util.function.ToDoubleFunction;

import gedi.util.io.text.jhp.TemplateGenerator;

public interface AutoTemplateGenerator<P,T extends AutoTemplateGenerator<P,T>> extends TemplateGenerator<T>, ToDoubleFunction<P> {

}
