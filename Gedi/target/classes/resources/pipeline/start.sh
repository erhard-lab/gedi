#!/bin/bash


<?JS for each (var n in js.getVariables(false).keySet()) {
var ff = new File(js.getVariable(n));
if (js.getVariable(n).class==String.class && ff.isAbsolute()) { ?>
export <? n ?>="<? print(js.getVariable(n)) ?>"
<?JS }} ?>
