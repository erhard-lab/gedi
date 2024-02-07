package gedi.grand3.estimation;

import java.io.IOException;
import java.util.HashMap;

import gedi.core.region.ImmutableReferenceGenomicRegion;
import gedi.util.FileUtils;
import gedi.util.datastructure.collections.intcollections.IntIterator;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;
import gedi.util.math.stat.inference.mixture.BiMixtureModelResult;

public class TargetEstimationResult implements BinarySerializable {

	
	public enum ModelType {
		Binom,TbBinom,TbBinomShape
	}
	
	private ImmutableReferenceGenomicRegion<String> target;
	private String category;
	
	private HashMap<Integer,Double> count = new HashMap<Integer, Double>();
	private HashMap<Integer,SingleEstimationResult>[] map;
	
	public TargetEstimationResult() {} // for deserialization
	
	public TargetEstimationResult(ImmutableReferenceGenomicRegion<String> target, String category, int numLabels) {
		this.target = target;
		this.category = category;
		map = new HashMap[numLabels];
		for (int i=0; i<numLabels; i++)
			map[i] = new HashMap<Integer, TargetEstimationResult.SingleEstimationResult>();
	}

	public ImmutableReferenceGenomicRegion<String> getTarget() {
		return target;
	}
	public String getCategory() {
		return category;
	}
	public String getGenome() {
		if (category.startsWith("Exonic")) return category.substring(8,category.length()-1);
		if (category.startsWith("Intronic")) return category.substring(10,category.length()-1);
		throw new RuntimeException("Invalid category format!");
	}
	public boolean isIntronic() {
		return category.startsWith("Intronic");
	}
	public boolean isExonic() {
		return category.startsWith("Exonic");
	}

	public void set(int i, double count) {
		this.count.put(i,count);
	}
	
	public SingleEstimationResult get(int t, int i) {
		return map[t].get(i);
	}
	
	public IntIterator iterateCounts() {
		return EI.wrap(count.keySet()).castInt();
	}
	
	public IntIterator iterate(int t) {
		return EI.wrap(map[t].keySet()).castInt();
	}
	
	public double getCountOrZero(int i) {
		Double r = count.get(i);
		if (r==null) return 0;
		return r;
	}
	
	public int getNumNonzeroCounts() {
		return count.size();
	}
	
	public int getNumNonzeroNtrs(int t) {
		return map[t].size();
	}
	
	
	public void setTargetShapeInfo(int t, int i, double shape, double logLikShape, double logLikGlobal) {
		SingleEstimationResult est = map[t].get(i);
		if (est==null) map[t].put(i, est = new SingleEstimationResult());
		est.shape = new ShapeEstimationResult(shape, logLikShape, logLikGlobal);
	}
	
	public void setTargetEstimateBinom(int t, int i, BiMixtureModelResult result) {
		SingleEstimationResult est = map[t].get(i);
		if (est==null) map[t].put(i, est = new SingleEstimationResult());
		est.binomMix = result;
	}
	public void setTargetEstimateTbBinom(int t, int i, BiMixtureModelResult result) {
		SingleEstimationResult est = map[t].get(i);
		if (est==null) map[t].put(i, est = new SingleEstimationResult());
		est.tbbinomMix = result;
	}
	public void setTargetEstimateTbBinomShape(int t, int i, BiMixtureModelResult result) {
		SingleEstimationResult est = map[t].get(i);
		if (est==null) map[t].put(i, est = new SingleEstimationResult());
		est.tbbinomShapeMix = result;
	}
	
	
	public boolean isEmpty() {
		return EI.wrap(map).mapToInt(m->m.size()).sum()==0;
	}


	@Override
	public void serialize(BinaryWriter out) throws IOException {
		FileUtils.writeReferenceSequence(out, target.getReference());
		FileUtils.writeGenomicRegion(out, target.getRegion());
		out.putString(target.getData());
		out.putString(category);
		
		out.putCInt(count.size());
		for (Integer c : count.keySet()) {
			out.putCInt(c);
			out.putDouble(count.get(c));
		}
		
		out.putCInt(map.length);

		for (int t=0; t<map.length; t++) {
			out.putCInt(map[t].size());
			for (Integer i : map[t].keySet()) {
				out.putCInt(i);
				map[t].get(i).serialize(out);
			}
		}
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		target = new ImmutableReferenceGenomicRegion<>(
				FileUtils.readReferenceSequence(in),
				FileUtils.readGenomicRegion(in),
				in.getString()
				);
		category = in.getString();
		
		int numCond = in.getCInt();
		for (int i=0; i<numCond; i++) 
			count.put(in.getCInt(), in.getDouble());
		
		int numLabels = in.getCInt();
		map = new HashMap[numLabels];
		for (int i=0; i<numLabels; i++)
			map[i] = new HashMap<Integer, TargetEstimationResult.SingleEstimationResult>();
		
		
		for (int t=0; t<map.length; t++) {
			int size = in.getCInt();
		
			for (int i=0; i<size; i++) {
				SingleEstimationResult s = new SingleEstimationResult();
				map[t].put(in.getCInt(), s);
				s.deserialize(in);
			}
		}
	}

	public static class SingleEstimationResult implements BinarySerializable{
		private ShapeEstimationResult shape;
		private BiMixtureModelResult binomMix;
		private BiMixtureModelResult tbbinomMix;
		private BiMixtureModelResult tbbinomShapeMix;
		
		public SingleEstimationResult() {}

		public ShapeEstimationResult getShape() {
			return shape;
		}
		
		public BiMixtureModelResult getModel(ModelType type) {
			switch (type) {
			case Binom: return binomMix;
			case TbBinom: return tbbinomMix;
			case TbBinomShape: return tbbinomShapeMix;
			default: throw new RuntimeException("Unknown model type!");
			}
		}

		public BiMixtureModelResult getBinomMix() {
			return binomMix;
		}

		public BiMixtureModelResult getTbbinomMix() {
			return tbbinomMix;
		}

		public BiMixtureModelResult getTbbinomShapeMix() {
			return tbbinomShapeMix;
		}

		@Override
		public void serialize(BinaryWriter out) throws IOException {
			int flag = 0;
			if (shape!=null) flag|=1<<0;
			if (binomMix!=null) flag|=1<<1;
			if (tbbinomMix!=null) flag|=1<<2;
			if (tbbinomShapeMix!=null) flag|=1<<3;
			if (flag!=15)
				throw new RuntimeException();
			out.putByte(flag);
			
			if (shape!=null) shape.serialize(out);
			if (binomMix!=null) binomMix.serialize(out);
			if (tbbinomMix!=null) tbbinomMix.serialize(out);
			if (tbbinomShapeMix!=null) tbbinomShapeMix.serialize(out);
			
		}

		@Override
		public void deserialize(BinaryReader in) throws IOException {
			int flags = in.getByte();
			if (flags!=15)
				throw new RuntimeException();
			if ((flags & (1<<0))!=0) {
				shape = new ShapeEstimationResult();
				shape.deserialize(in);
			}
			
			if ((flags & (1<<1))!=0) {
				binomMix = new BiMixtureModelResult();
				binomMix.deserialize(in);
			}
			
			if ((flags & (1<<2))!=0) {
				tbbinomMix = new BiMixtureModelResult();
				tbbinomMix.deserialize(in);
			}
			
			if ((flags & (1<<3))!=0) {
				tbbinomShapeMix = new BiMixtureModelResult();
				tbbinomShapeMix.deserialize(in);
			}
		}
		
		
		
	}

	public static class ShapeEstimationResult implements BinarySerializable {
		private double shape;
		private double logLikShape;
		private double logLikGlobal;
		
		public ShapeEstimationResult() {}
		public ShapeEstimationResult(double shape, double logLikShape, double logLikGlobal) {
			this.shape = shape;
			this.logLikShape = logLikShape;
			this.logLikGlobal = logLikGlobal;
		}
		@Override
		public String toString() {
			return "shape=" + shape + ", logLikShape=" + logLikShape + ", logLikGlobal="
					+ logLikGlobal;
		}
		public double getShape() {
			return shape;
		}
		public double getLogLikShape() {
			return logLikShape;
		}
		public double getLogLikGlobal() {
			return logLikGlobal;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			long temp;
			temp = Double.doubleToLongBits(logLikGlobal);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(logLikShape);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			temp = Double.doubleToLongBits(shape);
			result = prime * result + (int) (temp ^ (temp >>> 32));
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			ShapeEstimationResult other = (ShapeEstimationResult) obj;
			if (Double.doubleToLongBits(logLikGlobal) != Double.doubleToLongBits(other.logLikGlobal))
				return false;
			if (Double.doubleToLongBits(logLikShape) != Double.doubleToLongBits(other.logLikShape))
				return false;
			if (Double.doubleToLongBits(shape) != Double.doubleToLongBits(other.shape))
				return false;
			return true;
		}
		@Override
		public void serialize(BinaryWriter out) throws IOException {
			out.putDouble(shape);
			out.putDouble(logLikShape);
			out.putDouble(logLikGlobal);
		}
		@Override
		public void deserialize(BinaryReader in) throws IOException {
			shape=in.getDouble();
			logLikShape=in.getDouble();
			logLikGlobal=in.getDouble();
		}
		
		
	}

	

	

	
	
}
