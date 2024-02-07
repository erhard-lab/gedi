package gedi.util.math.stat.kernel;

import java.util.function.DoubleBinaryOperator;
import java.util.function.DoubleUnaryOperator;

public interface Kernel extends DoubleUnaryOperator {
	double halfSize();
	default String name() {
		String n = getClass().getSimpleName();
		if (n.endsWith("Kernel")) n = n.substring(0, n.length()-"Kernel".length());
		return n;
	}
	
	default PreparedIntKernel prepare() {
		return new PreparedIntKernel(this,false);
	}
	
	default PreparedIntKernel prepare(boolean normalize) {
		return new PreparedIntKernel(this,true);
	}
	
	/**
	 * first operand is the kernel operand, second is the result of the parent kernel
	 * @param parent
	 * @param op
	 */
	default Kernel decorate(DoubleBinaryOperator op) {
		return new DecoratedKernel(this, op);
	}
	
	
}
