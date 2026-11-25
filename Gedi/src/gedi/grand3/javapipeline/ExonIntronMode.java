package gedi.grand3.javapipeline;

public enum ExonIntronMode {
	Standard{
		public boolean getExonForGlobal() {return true;	}
		public boolean getIntronForGlobal() {return false;}
		public boolean getExonForTarget() {return true;}
		public boolean getIntronForTarget() {return false;}
	},
	Talkseq {
		public boolean getExonForGlobal() {return false;	}
		public boolean getIntronForGlobal() {return true;}
		public boolean getExonForTarget() {return true;}
		public boolean getIntronForTarget() {return true;}
	}, 
	TalkseqEE {
		public boolean getExonForGlobal() {return true;	}
		public boolean getIntronForGlobal() {return false;}
		public boolean getExonForTarget() {return true;}
		public boolean getIntronForTarget() {return false;}
	}, 
	TalkseqIE {
		public boolean getExonForGlobal() {return false;	}
		public boolean getIntronForGlobal() {return true;}
		public boolean getExonForTarget() {return true;}
		public boolean getIntronForTarget() {return false;}
	}, 
	Intronic {
		public boolean getExonForGlobal() {return false;	}
		public boolean getIntronForGlobal() {return true;}
		public boolean getExonForTarget() {return false;}
		public boolean getIntronForTarget() {return true;}
	};
	
	public abstract boolean getExonForGlobal();
	public abstract boolean getIntronForGlobal();
	public abstract boolean getExonForTarget();
	public abstract boolean getIntronForTarget();
	
}