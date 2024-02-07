package gedi.grand3.targets;

public class CompatibilityCategory {
	
	private String name;
	private int id;
	private boolean inferClipping;
	private boolean estimateGlobal;
	private boolean estimateTarget;
	private boolean reverseLabel;
	

	public CompatibilityCategory(String name, int id, boolean inferClipping, boolean estimateGlobal,
			boolean estimateTarget, boolean reverseLabel) {
		this.name = name;
		this.id = id;
		this.inferClipping = inferClipping;
		this.estimateGlobal = estimateGlobal;
		this.estimateTarget = estimateTarget;
		this.reverseLabel = reverseLabel;
	}

	public String getName() { return name; }
	public int id() { return id; }
	public boolean useToInferClipping() { return inferClipping; }
	public boolean useToEstimateGlobalParameters() { return estimateGlobal; }
	public boolean useToEstimateTargetParameters() { return estimateTarget; }
	public boolean reversesLabel()  { return reverseLabel; }

	
	@Override
	public String toString() {
		return name;
	}
}