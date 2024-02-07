package gedi.gui.gtracks.style;

import java.util.function.Function;

import gedi.util.PaintUtils;
import gedi.util.dynamic.DynamicObject;

public interface GTrackStyleParser extends Function<DynamicObject,Object> {
	
	
	public static class GTrackNameStyleParser implements GTrackStyleParser {
		@Override
		public Object apply(DynamicObject t) {
			return t.asString();
		}
	}

	
	public static class GTrackColorStyleParser implements GTrackStyleParser {
		@Override
		public Object apply(DynamicObject t) {
			return PaintUtils.parseColor(t.asString());
		}
	}

	public static class GTrackFillStyleParser implements GTrackStyleParser {
		@Override
		public Object apply(DynamicObject t) {
			return PaintUtils.parseColor(t.asString());
		}
	}
	
	public static class GTrackSizeStyleParser implements GTrackStyleParser {
		@Override
		public Object apply(DynamicObject t) {
			return t.asDouble();
		}
	}

	
}