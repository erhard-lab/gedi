package gedi.centeredDiskIntervalTree;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.function.Supplier;

import gedi.core.region.GenomicRegion;
import gedi.util.FileUtils;
import gedi.util.dynamic.DynamicObject;
import gedi.util.io.randomaccess.PageFile;
import gedi.util.io.randomaccess.PageFileWriter;
import gedi.util.orm.BinaryBlob;

public class InternalInMemoryCenteredDiskIntervalTreeBuilder<D> extends CenteredDiskIntervalTreeBuilder<D> {

	public static final String MAGIC = "CITI";
	
	private BinaryBlob data = new BinaryBlob();
	
	private long offset = 0;
	
	public InternalInMemoryCenteredDiskIntervalTreeBuilder(String prefix, DynamicObject globalInfo) throws IOException {
		this(System.getProperty("java.io.tmpdir"),prefix, globalInfo);
	}
	
	
	public InternalInMemoryCenteredDiskIntervalTreeBuilder(String tmpFolder, String prefix, DynamicObject globalInfo) throws IOException {
		super(true,MAGIC,prefix,tmpFolder);
	}
	
	@Override
	public void toDisk() throws IOException {
		
	}
	
	public void add(GenomicRegion region, D data) throws IOException {
		long ptr = this.data.position()+offset;
		if (region.getBoundary(0)>Integer.MAX_VALUE/2 || region.getBoundary(0)<0)
			return;
		
		this.data.putCInt(region.getNumParts());
		int start = region.getBoundary(0);
		this.data.putCInt(start);
		for (int i=1; i<region.getNumBoundaries(); i++)
			this.data.putCInt(region.getBoundary(i)-start);
		
		FileUtils.serialize(data,this.data);
		
//		data.serialize(this.data);
		add(region,ptr);
	}
	
	public InternalInMemoryCenteredDiskIntervalTreeBuilder<D> build(PageFileWriter out) throws IOException {
		super.build(out);
		
		data.finish(false);
		while (!data.eof()) {
			out.put(data.get());
		}
		data.close();
		
		return this;
	}
	
	
}
