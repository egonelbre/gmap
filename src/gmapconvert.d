module gmapconvert;

import std.stdio;
import std.file;
import std.math;
import std.stream;
import std.string;

import gmap.gmap;

import tools.flags;

void help(int idx, string value) {
    writef(
"Usage: gmapconvert -i <input> <dataset>

Program converts .ped and .map file to gmap format.

Input:
    There must be files <input>.ped, <input>.map,
    <input>.gt

Additional Info:
    <input>.gt must contain all genotypes used in <input>.ped.
    Each line in <input>.gt specifies a genotype.
    If some genotypes can be considered equivalent (such
    as AB and BA) they can be on the same line (tab delimited).
    Afterwards both are considered equivalent and can be
    used interexchangably.

Output:
    <dataset>.gmap : binary genotype data
    <dataset>.gidx : genotypes index
    <dataset>.midx : marker index
    <dataset>.pidx : people index
");
    abort_program();
};

string prefix = "";
string prefixOut = "";

int main(string[] args)
{
    flagsNoArg["-h"]  = &help;
    flags["-i"] =
        function void(int idx, string value){ prefix = value; };
    flags["def"] =
        function void(int idx, string value){
            if( prefixOut == "")
                prefixOut = value;
            else{
                writefln("Invalid argument %s", value);
                help(0,"");
            }
        };
    parseFlags(args);

    if(prefix == "" || prefixOut == ""){
        writefln("No file specified.");
        help(0,"");
    }
    
    string pedfile = prefix ~ ".ped";
    string mapfile = prefix ~ ".map";
    string genotypefile = prefix ~ ".gt";
    
    // input files
    auto ped = new BufferedFile(pedfile);
    auto map = new BufferedFile(mapfile);
    auto gt  = new BufferedFile(genotypefile);
    
    // our output files
    auto gmap  = new BufferedFile(prefixOut ~ sufGMap, FileMode.OutNew);
    auto pidx  = new BufferedFile(prefixOut ~ sufPeople, FileMode.OutNew);
    auto midx  = new BufferedFile(prefixOut ~ sufMarkers, FileMode.OutNew);
    auto gidx  = new BufferedFile(prefixOut ~ sufGenotypes, FileMode.OutNew);
    
    // create to ".gidx" mapping for genotype -> idx
    writefln("creating .gidx ...");
    int[char[]] genotypemap;
    int genotypecount = 0;
    while( !gt.eof() ){
        string genotype = gt.readLine();
        string[] genotypes = split(genotype, "\t");
        for(int i = 0; i < genotypes.length; i++)
            genotypemap[genotypes[i]] = genotypecount;
        genotypecount++;
    }
    gidx.writefln(genotypecount);
    foreach(genotype, index; genotypemap)
        gidx.writefln("%s\t%d", genotype, index);
    
    gidx.close();
    
    writefln("creating .pidx ...");
    ulong[] plate = new ulong[2100]; // will be resized if nessecary
    string[] platename = new string[2100];
    ulong iplate = 0;
    
    while (!ped.eof()){       
        //the syntax of the file
        //1#Geenivaramu_plaat1#A01	1#Geenivaramu_plaat1#A01	0	0	0	0	...\n
        //                                  this position will be remembered----^
        char c;
        string buf;
        
        for(int i=0; i<6; i++){
            ped.read(c);
            while( c != '\t' ){
                buf ~= c; 
                ped.read(c);
            }
            if (i < 5)
                buf ~= c;
        }
        
        plate[iplate] = ped.position;
        platename[iplate] = buf;
        iplate++;
        
        if (plate.length <= iplate)
            plate.length = plate.length << 2;
        if (platename.length <= iplate)
            platename.length = platename.length << 2;
        
        while( c != '\n' )
            ped.read(c);
    }
    
    plate.length = iplate;
    platename.length = iplate;
    
    pidx.writefln(platename.length);
    foreach(index, name; platename)
        pidx.writefln("%d\t%s", index, name);
    pidx.close();
    
    // if we need appending plates, we need to make a empty area to 
    // the end
    ulong count = plate.length / intgbits;
    if (plate.length % intgbits != 0)
        count++;

    // for speed we may need to construct multiple linebuffers
    intg[][] linebuffer   = new intg[][genotypecount];
    for(int i = 0; i < genotypecount; i++)
        linebuffer[i] = new intg[count];
    
    // this is how many markers will be rotated simultaniously
    ulong pedbuffercount = 100;
    
    byte[][] pedbuffer = new byte[][pedbuffercount];
    for(int i = 0; i < pedbuffercount; i++)
        pedbuffer[i] = new byte[plate.length];
    
    writefln("converting data...");
    
    string marker_format(string inp){
        string[] inps = split(inp);
        return format("%s\t%s\t%s\t%s", inps[0], inps[1], inps[2], inps[3]);
    }
    
    // count then number of markers
    // assumes each line is a marker definition
    ulong marker_count = 0;
    foreach(string line; map){
        if ( line != "" )
            marker_count++;
    }
    midx.writefln("%d", marker_count);
    
    map.position = 0;
    string marker = marker_format(map.readLine());
    
    CONVERT_LOOP:
    while( marker != "" ){
        debug writefln("reading marker %s", marker);
        {   // fill up our pedbuffer, this rotates the marker column
            char c;
            char[100] buf;
            
            for(ulong i = 0; i < plate.length; i++){
                if (i == 0)
                    ped.seek(plate[i], SeekPos.Set);
                else
                    ped.seek(plate[i] - ped.position, SeekPos.Current);
                    
                for(ulong u = 0; u < pedbuffercount; u++){
                    ulong k = 0;
                    ped.read(c);
                    while( c != '\t' && c != '\n'){
                        buf[k] = c; k++;
                        if (k >= 100) k = 0;
                        ped.read(c);
                    }
                    scope(failure){ 
                        writefln("Error accessing: '%s' N'%s' @ %d - [%d][%d].%d", 
                                    buf[0..k], 
                                    c, 
                                    ped.position,
                                    u, i, plate[i]);
                    }
                    pedbuffer[u][i] = genotypemap[buf[0..k]];
                    if (c == '\n')
                        break;
                }
                
                plate[i] = ped.position;
            }
        }
              
        for(ulong u = 0; u < pedbuffercount; u++){
            ulong p = 0;  // position in intg
            ulong k = 0;  // position in singlebuffer
            
            intg[] singlebuffer = new intg[genotypecount];
            
            // this func copies singlebuffer to linebuffer
            void flush_singlebuffer(){
                for(ulong n = 0; n < genotypecount; n++ )
                        linebuffer[n][k] = singlebuffer[n];
            }            
            
            for(ulong i = 0; i < pedbuffer[u].length; i++){
                byte n = pedbuffer[u][i];
                singlebuffer[n] = singlebuffer[n] | (1L << p);
                
                if (p == intgbits - 1){
                    flush_singlebuffer();
                    k++;
                    singlebuffer[] = 0;
                }
                p = (p + 1) % intgbits;
            }
            // if the last singlebuffer wasn't flushed
            if (p <> 0)
                flush_singlebuffer();
            
            midx.writef("%s", marker);
            midx.writef("\t%d", gmap.position);
            for(ulong n = 0; n < genotypecount; n++ )
                gmap.writeBlock(linebuffer[n].ptr, intg.sizeof * count);
                
            midx.writefln();
            
            if( map.eof() )
                break CONVERT_LOOP;
            marker = marker_format(map.readLine());
        }
    }
    gmap.close();
    midx.close();
    
    return 1;
}
