package gedi.util.oml.petrinet;

import gedi.core.data.mapper.GenomicRegionDataMapper;
import gedi.core.data.mapper.GenomicRegionDataMappingJob;
import gedi.core.data.mapper.MutableDemultiplexMapper;
import gedi.core.genomic.Genomic;
import gedi.gui.genovis.VisualizationTrack;
import gedi.util.ReflectionUtils;
import gedi.util.StringUtils;
import gedi.util.functions.EI;
import gedi.util.io.text.LineIterator;
import gedi.util.job.PetriNet;
import gedi.util.job.Place;
import gedi.util.job.Transition;
import gedi.util.mutable.Mutable;
import gedi.util.oml.OmlInterceptor;
import gedi.util.oml.OmlNode;
import gedi.util.oml.OmlNodeExecutor;
import gedi.util.oml.OmlReader;
import gedi.util.oml.cps.CpsReader;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.xml.sax.InputSource;

import executables.Template;

@SuppressWarnings({"unchecked","rawtypes"})
public class Pipeline implements OmlInterceptor {

	
	private static final String INPUT_ATTRIBUTE = "input";
	private static final String INPUT_SETTER = "setInput";
	
	private ArrayList<VisualizationTrack> tracks = new ArrayList<VisualizationTrack>();
	private PetriNet pn;
	
//	private Place lastPlace;
	private HashMap<String, Place> idToPlace = new HashMap<String, Place>();

	private Place[] inputs;
	private String inputIds;
	private Genomic genomic;
	
	public Pipeline() {
		pn = new PetriNet();
	}

	public PetriNet getPetriNet() {
		if (!pn.isPrepared())
			pn.prepare();
		return pn;
	}

	public Genomic getGenomic() {
		return genomic;
	}
		
	
	public void sortPlusMinusTracks() {
		Collections.sort(tracks,(a,b)->pmPref(a).compareTo(pmPref(b)));
	}
	
	private static String pmPref(VisualizationTrack t) {
		String p = t.getId();
//		if (p.contains("4sU.Mock") || p.contains("4sU.8hpi")) {
//			if (p.startsWith("+")) return "-";
//			if (p.startsWith("-")) return "+";
//		}
		if (p.startsWith("+")) return "+";
		if (p.startsWith("-")) return "-";
		return "";
	}
	
	
	public ArrayList<VisualizationTrack> getTracks() {
		return tracks;
	}
	
	@Override
	public boolean useForSubtree() {
		return true;
	}
	
	private HashMap<Class,ArrayList> objects = new HashMap<Class, ArrayList>();
	
	public <T> ArrayList<T> getObjects(Class<T> cls) {
		if (!objects.containsKey(cls)) return new ArrayList();
		return objects.get(cls);
	}

	public void add(GenomicRegionDataMapper o) {
		if (o instanceof VisualizationTrack)
			tracks.add((VisualizationTrack<?,?>) o);
	}
	
	@Override
	public void setObject(OmlNode node, Object obj, String id, String[] classes,
			HashMap<String, Object> context) {
		
		if (obj instanceof Genomic && genomic==null)
			genomic = (Genomic) obj;
		
		if (obj instanceof GenomicRegionDataMapper) {
			
			GenomicRegionDataMapper o = (GenomicRegionDataMapper) obj;
			
			ArrayList l = objects.get(o.getClass());
			if (l==null) objects.put(o.getClass(), l = new ArrayList());
			l.add(o);
			
			GenomicRegionDataMappingJob job = new GenomicRegionDataMappingJob(o);
			
			if (job.isTupleInput()) 
				job.setTupleSize(inputs.length);
			
			Transition t = pn.createTransition(job);
			Place p = pn.createPlace(job.getOutputClass());
			pn.connect(t, p);
			
			idToPlace.put(id,p);
			((GenomicRegionDataMappingJob) p.getProducer().getJob()).setId(id);
			
			if (inputs!=null) {
				for (int i=0; i<inputs.length; i++)
					pn.connect(inputs[i], t, i);
			}
			if (inputIds!=null) {
				Method setter = ReflectionUtils.findMethod(o, INPUT_SETTER,String.class);
				if (setter!=null)
					try {
						setter.invoke(o, inputIds);
					} catch (IllegalAccessException | IllegalArgumentException
							| InvocationTargetException e) {
						throw new RuntimeException("Could not execute setter: "+setter,e);
					}
			}
		}
	}
	
	
	
	@Override
	public LinkedHashMap<String, String> getAttributes(OmlNode node, LinkedHashMap<String, String> attributes,
			HashMap<String, Object> context) {
		
		if (context.containsKey(OmlNodeExecutor.INLINED_CALL)) return attributes;
		
		inputs = null;
		inputIds = null;
		
		LinkedHashMap<String, String> re = attributes;
		if (re.containsKey(INPUT_ATTRIBUTE)) {
			re = new LinkedHashMap<String, String>(re);
			this.inputIds = re.get(INPUT_ATTRIBUTE);
			String[] inputIds = StringUtils.split(this.inputIds,',');
			inputs = new Place[inputIds.length];
			for (int i=0; i<inputs.length; i++) {
				inputs[i] = idToPlace.get(inputIds[i]);
				int mi = getMutableIndex(inputIds[i]);
				if (inputs[i]==null && mi!=-1) {
					
					// add <MutableDemultiplexMapper input="pref" index="mi" id="pref[mi]" />
					
					String pref = inputIds[i].substring(0, inputIds[i].lastIndexOf('['));
					String full = inputIds[i];
					if (idToPlace.containsKey(pref) && Mutable.class.isAssignableFrom(idToPlace.get(pref).getTokenClass())) {
						MutableDemultiplexMapper m = new MutableDemultiplexMapper(mi);
						
						Place[] inputssave = this.inputs;
						String inputIdssave = this.inputIds;
						this.inputIds = pref;
						this.inputs = new Place[] {idToPlace.get(pref)};
						
						setObject(null, m, full, null, null);
						
						this.inputIds = inputIdssave;
						this.inputs = inputssave;
						this.inputs[i] = idToPlace.get(inputIds[i]);
						
					}
				}
				if (inputs[i]==null)
					throw new RuntimeException("Input with id "+inputIds[i]+" unknown!");
			}
			
			
			
			re.remove(INPUT_ATTRIBUTE);
		}
		
		return re;
	}

	private int getMutableIndex(String id) {
		if (!id.endsWith("]")) return -1;
		int s = id.lastIndexOf('[');
		if (s==-1) return -1;
		String n = id.substring(s+1, id.length()-1);
		if (StringUtils.isInt(n))
			return Integer.parseInt(n);
		return -1;
	}

	
	public static Pipeline fromOml(String fileResourceOrSource) throws IOException {
		String cps = new LineIterator(Pipeline.class.getResourceAsStream("/resources/colors.cps")).concat("\n");
		OmlNodeExecutor oml = new OmlNodeExecutor()
				.addInterceptor(new CpsReader().parse(cps));
		
		if (new File(fileResourceOrSource).exists())
			return (Pipeline)oml.execute(new OmlReader().parse(new File(fileResourceOrSource)));
		InputStream res = Pipeline.class.getResourceAsStream(fileResourceOrSource);
		if (res!=null)
			return (Pipeline)oml.execute(new OmlReader().parse(new InputSource(res)));
		
		return (Pipeline)oml.execute(new OmlReader().parse(fileResourceOrSource));
	}

		
}
