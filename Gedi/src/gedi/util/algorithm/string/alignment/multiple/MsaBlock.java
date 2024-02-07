package gedi.util.algorithm.string.alignment.multiple;

import java.io.IOException;
import java.util.LinkedHashMap;

import gedi.app.extension.GlobalInfoProvider;
import gedi.util.dynamic.DynamicObject;
import gedi.util.functions.EI;
import gedi.util.io.randomaccess.BinaryReader;
import gedi.util.io.randomaccess.BinaryWriter;
import gedi.util.io.randomaccess.serialization.BinarySerializable;

public class MsaBlock implements BinarySerializable, GlobalInfoProvider {

	public static final String ROWSATTRIBUTE = "ROWCOUNT";
	
	private char[][] rows;
	
	public MsaBlock(){}
	public MsaBlock(char[][] rows) {
		this.rows = rows;
	}
	
	public char[][] getRows() {
		return rows;
	}

	@Override
	public void serialize(BinaryWriter out) throws IOException {
		DynamicObject gi = out.getContext().getGlobalInfo();
		if (!gi.hasProperty(ROWSATTRIBUTE))
			out.putCInt(rows.length);
		out.putCInt(rows[0].length);
		
		
		for (int c=0; c<rows[0].length; c++)
			for (int r=0; r<rows.length; r++)
				out.putAsciiChar(rows[r][c]);
		
	}

	@Override
	public void deserialize(BinaryReader in) throws IOException {
		int rows;
		
		DynamicObject gi = in.getContext().getGlobalInfo();
		if (!gi.hasProperty(ROWSATTRIBUTE))
			rows = in.getCInt();//conditions
		else
			rows = gi.getEntry(ROWSATTRIBUTE).asInt();
		
		int cols = in.getCInt();
		
		this.rows = new char[rows][cols];
		for (int c=0; c<cols; c++) 
			for (int r=0; r<rows; r++)
				this.rows[r][c] = in.getAsciiChar();
		
	}
	
	@Override
	public String toString() {
		return EI.wrap(rows).map(c->String.valueOf(c)).concat("\n");
	}
	@Override
	public DynamicObject getGlobalInfo() {
		LinkedHashMap<String, Object> map = new LinkedHashMap<String, Object>();
		map.put(ROWSATTRIBUTE, rows.length);
		return DynamicObject.from(map);
	}
	
}
