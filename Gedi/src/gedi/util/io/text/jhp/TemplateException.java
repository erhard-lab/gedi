package gedi.util.io.text.jhp;

public class TemplateException extends RuntimeException {

	public TemplateException(String msg) {
		super(msg);
	}
	
	public TemplateException(String msg, Throwable e) {
		super(msg,e);
	}
	
}
