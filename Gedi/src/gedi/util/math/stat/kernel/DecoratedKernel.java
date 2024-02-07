package gedi.util.math.stat.kernel;

import java.util.function.DoubleBinaryOperator;

public class DecoratedKernel implements Kernel {
	
	private Kernel parent;
	private DoubleBinaryOperator op;
	
	/**
	 * first operand is the kernel operand, second is the result of the parent kernel
	 * @param parent
	 * @param op
	 */
	public DecoratedKernel(Kernel parent, DoubleBinaryOperator op) {
		this.parent = parent;
		this.op = op;
	}

	@Override
	public double applyAsDouble(double operand) {
		return op.applyAsDouble(operand,parent.applyAsDouble(operand));
	}

	@Override
	public double halfSize() {
		return parent.halfSize();
	}

}
