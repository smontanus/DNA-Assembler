/* 
* dnaassembler.js
*
* Functions for DNA Assembler app. 
*
*/
//Waiting on response spinner.
function loading(fqfile){
    $(".loading").show();
    $(".block").hide();
    window.location = fqfile;
}