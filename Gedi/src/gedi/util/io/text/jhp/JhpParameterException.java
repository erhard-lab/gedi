package gedi.util.io.text.jhp;


public class JhpParameterException extends RuntimeException {

	private static final long serialVersionUID = 227393818068359180L;
	private String variableName;

	
	public JhpParameterException(String variableName) {
		this(variableName,null);
	}

	public JhpParameterException(String variableName,Throwable cause) {
		super("Unknown variable: "+variableName,cause);
		this.variableName = variableName;
	}

	public String getVariableName() {
		return variableName;
	}
	
	
}
