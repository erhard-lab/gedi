package gedi.util.parsing;

import gedi.core.data.reads.ReadCountMode;
import gedi.util.PaintUtils;

import java.awt.Color;
import java.awt.Paint;

public class ReadCountModeParser implements Parser<ReadCountMode> {

	@Override
	public ReadCountMode apply(String t) {
		return ReadCountMode.valueOf(t);
	}

	@Override
	public Class<ReadCountMode> getParsedType() {
		return ReadCountMode.class;
	}

}
