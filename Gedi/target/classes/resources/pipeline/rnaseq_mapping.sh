<?JS0

varin("wd","Working directory",true);
varin("tmp","Temp directory",true);
varin("tokens","Array of pipeline tokens (to resolve dependencies of programs)",true);
varin("name","Name for output files",true);
varin("datasets","Dataset definitions",true);
varin("references","Definitions of reference sequences",true);
varin("test","Test with the first 10k sequences",false);
varin("maxpar","Maximum number of parallel threads",false);
varin("root","root folder for all fastqs (fastq or fastq-regex)",false);
varin("keepnoumi","Also keep the cit file without umi deduplication",false);
varin("sharedMem","Load STAR index into shared mem first!",false);

varout("reads","File name containing read mappings");

?>
<?JS
var getPrio=function(r) ParseUtils.parseEnumNameByPrefix(references[r], true, ReferenceType.class).prio;
var garray = EI.wrap(DynamicObject.from(references).getProperties()).filter(function(i) getPrio(i)>1).toArray();
var genomes = EI.wrap(garray).concat(" ");


var SLAMrescue;
var SLAMrescueIndex=SLAMrescue?ReadMapper.STAR.getIndex(Genomic.get(SLAMrescue),null):null
var starindex;
if (!starindex && garray.length==1) {
	starindex=ReadMapper.STAR.getIndex(Genomic.get(garray[0]),null);
	println("export starindex=\""+starindex+"\"");
	starindex="$starindex"
}	
if (!starindex) throw new RuntimeException("Specify a starindex! You can create one by gedi -e GenomicUtils -p -m star -g "+EI.wrap(garray).concat(" "));


// extract vars for genome indices
var discardGenomes = new LinkedHashMap();
for each (var index in EI.wrap(DynamicObject.from(references).getProperties()).filter(function(i) getPrio(i)==1).loop()) {
	var varname = "INDEX_"+index.replaceAll("[^A-z]","_");
	discardGenomes.put(index,"$"+varname);
	println("export "+varname+"=\""+ReadMapper.bowtie2.getIndex(Genomic.get(index),null)+"\"");
}


if (sharedMem && runner=="parallel") {?>


STAR --genomeLoad LoadAndExit --genomeDir <? starindex ?>
<?JS if (SLAMrescue) { ?>
STAR --genomeLoad LoadAndExit --genomeDir <? SLAMrescueIndex ?>
<?JS } 
}?>


<?JS
var sharedMem;
var root;
log.info("Root: "+root);

var resolve = function(f)  {
	if (root) {
		var pat = Pattern.compile("/"+f+"$").asPredicate();
		return EI.files(root,f.contains("/")).str().filter(pat);
	}
	return EI.singleton(new File(f).getAbsolutePath());
}
var maxpar;
if (!maxpar) {
	if (runner=="parallel") maxpar=6;
	else maxpar=1000000;
}

output.setExecutable(true);

var srrPat = Pattern.compile("SRR\\d+");
var nameToModeAndFiles = {};
var names = [];

var pairedend;
var adapter;
var trimmed;
var allfileslist=new ArrayList();

for each (var d in datasets) {
	if (d.hasOwnProperty("gsm")) {
		var allsrrs = new ArrayList();
		var allgsm = d.gsm;
		if ( typeof allgsm === 'string' ) allgsm = [allgsm];
		for each (var gsm in allgsm) {
			var jsson = new LineIterator(new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term="+gsm+"&retmode=json").openStream())
					.concat("");
			var json2 = DynamicObject.parseJson(jsson);
			var arr = json2.get(".esearchresult.idlist").asArray();
			if (arr.length!=1) throw new RuntimeException("Did not get a unique id for "+d.gsm+": "+json2);
			var xml = new LineIterator(new URL("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="+arr[0].asString()+"&rettype=docsum").openStream())
				.concat("");
			var srrs = EI.wrap(srrPat.matcher(xml)).sort().toArray(String.class);
			allsrrs.addAll(Arrays.asList(srrs));
			Thread.sleep(1000);
		}
		allsrrs = EI.wrap(allsrrs).sort().toArray(String.class);
		log.info("SRA entry: Name: "+d.name+" gsm: "+d.gsm+" id: "+arr[0].asString()+" SRR: "+Arrays.toString(allsrrs));
		nameToModeAndFiles[d.name] = ["SRR",allsrrs,d.hasOwnProperty("adapter")?d.adapter:adapter,false];
	} else if (d.hasOwnProperty("sra")) {
		var srrs = d.sra;
		log.info("SRA entry: Name: "+d.name+" srr: "+srrs);
		nameToModeAndFiles[d.name] = ["SRR",JS.array(srrs),d.hasOwnProperty("adapter")?d.adapter:adapter,false];
	} else if (d.hasOwnProperty("fastq")) {
		var fastq = d.fastq;
		var files = EI.wrap(JS.array(fastq)).unfold(resolve).toArray(String.class);
		for each (var file in files) {
			allfileslist.add(file);
			if (pairedend) allfileslist.add(FileUtils.findPartnerFile(file,pairedend));
		}
		log.info("Fastq entry: Name: "+d.name+" fastq: "+EI.wrap(files).concat(","));
		var missing = EI.wrap(files).filter(function(f) !new File(f).exists()).concat(",");
		if (missing.length()>0) throw new RuntimeException("Files are missing: <"+missing+"> (n="+missing.length()+")");
		if (files.length==0)  throw new RuntimeException("Files not found: <"+fastq+"> in "+root);
		nameToModeAndFiles[d.name] = ["FASTQ",files,d.hasOwnProperty("adapter")?d.adapter:adapter,d.hasOwnProperty("trimmed")?d.trimmed:trimmed];
	} else {
		throw new RuntimeException("Each entry must either have a gsm,sra or fastq field!")
	}
	
	nameToModeAndFiles[d.name].push(d.barcodes);
	
	if (names.indexOf(d.name)!=-1) throw new Exception("Name "+d.name+" exists more than once!")
	names.push(d.name);
}

if (new HashSet(allfileslist).size()!=allfileslist.size()) throw new Exception("Some files are used in more than one sample!")


var id = name;
var tokens = tokens?tokens:[];

for (var name in nameToModeAndFiles) {
	var mode = nameToModeAndFiles[name][0];
	var files = EI.wrap(nameToModeAndFiles[name][1]).concat(" ");
	adapter = nameToModeAndFiles[name][2];
	trimmed = nameToModeAndFiles[name][3];
	barcodes = nameToModeAndFiles[name][4];
	processTemplate("rnaseq_mapping1.sh",output.file.getParent()+"/"+name+".bash");
	if (tokens.length>=maxpar) {
		prerunner(name,tokens);
		tokens=[];
	} else {
		prerunner(name); 
	}
	print("$wd/scripts/"+name+".bash"); tokens.push(postrunner(name)); println(""); 
		
}
?>

	
<?JS prerunner(id+".merge",tokens) ?>
gedi -e MergeCIT -c <?JS wd ?>/<?JS id ?>.cit <?JS print(EI.wrap(JS.array(names)).map(function(f) wd+"/"+f+".cit").concat(" ")); ?>
<?JS var end=[postrunner(id+".merge")] ?> 

<?JS if (sharedMem) { ?>

STAR --genomeLoad Remove --genomeDir <? starindex ?>
<?JS if (SLAMrescue) { ?>
STAR --genomeLoad Remove --genomeDir <? SLAMrescueIndex ?>
<?JS } ?>

<?JS } ?>


<?JS if (keepnoumi) { ?> 
<?JS prerunner(id+".mergenodedup",tokens) ?>
gedi -e MergeCIT -c <?JS wd ?>/nodedup/<?JS id ?>.cit <?JS print(EI.wrap(JS.array(names)).map(function(f) wd+"/nodedup/"+f+".cit").concat(" ")); ?>
<?JS var end=[end[0],postrunner(id+".mergenodedup")] ?> 
<?JS } ?>


<?JS

name = id;
tokens = end;
var reads = wd+"/"+id+".cit";

log.info("Total number of fastq files used: "+allfileslist.size());

?> 
