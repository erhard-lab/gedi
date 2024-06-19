<?JS0

varin("wd","Working directory",true);
varin("tmp","Temp directory",true);
varin("tokens","Array of pipeline tokens (to resolve dependencies of programs)",true);
varin("name","Name for output files",true);
varin("numi","Count nUMIs",true);
varin("references","Definitions of reference sequences",true);
varin("reads","Filename containing read mappings",true);

?>

<?JS
output.setExecutable(true);

var tokens;
var getPrio=function(r) ParseUtils.parseEnumNameByPrefix(references[r], true, ReferenceType.class).prio;
var garray = EI.wrap(DynamicObject.from(references).getProperties()).filter(function(i) getPrio(i)>1).toArray();
var rrnaarray = EI.wrap(DynamicObject.from(references).getProperties()).filter(function(i) getPrio(i)==1).toArray();
var genomes = EI.wrap(garray).concat(" ");
	

var numi;

?>

<?JS prerunner(id+".grandslam",tokens) ?>gedi -t <?JS tmp ?> -e Slam -reads <?JS reads ?> -genomic <?JS genomes ?> -prefix <?JS wd ?>/grandslam/<?JS name ?> -plot <? if (numi) print("-nUMI"); ?> -D -modelall <?JS var grandslam = postrunner(id+".grandslam") ?> 

<?JS prerunner(id+".grandslam_t15",tokens) ?>gedi -t <?JS tmp ?> -e Slam -trim5p 15 -reads <?JS reads ?> -genomic <?JS genomes ?> -prefix <?JS wd ?>/grandslam_t15/<?JS name ?> -plot <? if (numi) print("-nUMI"); ?> -D -modelall <?JS var grandslam_t15 = postrunner(id+".grandslam_t15") ?> 

<?JS prerunner(id+".grand3",tokens) ?>gedi -mem 64G -e Grand3 -auto -reads <?JS reads ?> -genomic <?JS genomes ?> -prefix <?JS wd ?>/grand3/<?JS name ?> -D -profile <?JS var grand3 = postrunner(id+".grand3") ?> 


<?JS

name = id;
tokens = [grandslam,grandslam_t15,grand3];

?> 
