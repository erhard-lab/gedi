package executables;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Stack;
import java.util.logging.Level;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import gedi.util.StringUtils;
import gedi.util.functions.EI;

public class PepXml2Csv {

	
	public static void main(String[] args) {
		if (args.length!=1 || !new File(args[0]).exists()) {
			System.err.println("PepXml2Csv <XYZ.pep.xml>");
			System.exit(1);
		}
		try {
			SAXParser saxParser = SAXParserFactory.newInstance().newSAXParser();
			PepHandler handler = new PepHandler();
			System.out.println(EI.wrap(handler.map.keySet()).concat(","));
			
			saxParser.parse(new File(args[0]), handler);
		} catch (Exception e) {
			throw new RuntimeException("Cannot parse input "+args[0],e);
		}
		
	}
	
	
	private static class PepHandler extends DefaultHandler {
		
		private LinkedHashMap<String,String> map = new LinkedHashMap<>();
		
		public PepHandler() {
			map.put("spectrum", "");
			map.put("start_scan", "");
			map.put("end_scan", "");
			map.put("precursor_neutral_mass", "");
			map.put("assumed_charge", "");
			map.put("index", "");
			map.put("hit_rank", "");
			map.put("peptide", "");
			map.put("massdiff", "");
			map.put("Denovo", "");
			map.put("PeaksScore", "");
		}
		
		@Override
		public void startElement(String uri, String localName, String qName,
				Attributes attributes) throws SAXException {

			if (qName.equals("spectrum_query") || qName.equals("search_hit")) {
				if (qName.equals("spectrum_query"))
					for (String k : map.keySet())
						map.put(k, "");
				
				for (int i=0; i<attributes.getLength(); i++) {
					String name = attributes.getQName(i);
					String val = attributes.getValue(i);
					if (name.equals("protein"))
						map.put("Denovo", val.startsWith("|denovo_")?"1":"0");
					if (map.containsKey(name))
						map.put(name, val);
				}
			}
			else if (qName.equals("alternative_protein")) {
				for (int i=0; i<attributes.getLength(); i++) {
					String name = attributes.getQName(i);
					String val = attributes.getValue(i);
					if (name.equals("protein"))
						map.compute("Denovo", (k,v)->val.startsWith("|denovo_")&&!v.equals("0")?"1":"0");
				}
			}
			else if (qName.equals("modification_info")) {
				for (int i=0; i<attributes.getLength(); i++) {
					String name = attributes.getQName(i);
					String val = attributes.getValue(i);
					if (name.equals("modified_peptide"))
						map.put("peptide", val);
				}
			}
			else if (qName.equals("search_score")) {
				if (EI.seq(0,attributes.getLength()).filter(i->attributes.getQName(i).equals("name")).count()>0)
					map.put("PeaksScore", EI.seq(0,attributes.getLength()).filter(i->attributes.getQName(i).equals("value")).map(i->attributes.getValue(i)).getUniqueResult(true, true));
			}

		}
		
		@Override
		public void endElement(String uri, String localName, String qName)
				throws SAXException {
			if (qName.equals("spectrum_query"))
				System.out.println(EI.wrap(map.values()).concat(","));
		}
		
	}
	
}
