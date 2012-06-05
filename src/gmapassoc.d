module gmapassoc;

import std.file;
import std.stream;
import std.string;

import std.stdio;
import std.math;

import gmap.gmap;
import gmap.mask;
import gmap.statistics;
import gmap.utils;

import tools.flags;

void help(int idx, string value) {
    writef(
"Usage: gmapassoc -g <group.pset> <dataset>

Program does a basic case/control association analysis.

Input:
    There must be files <dataset>.gmap, <dataset>.gidx,
    <dataset>.midx, <dataset>.aidx and <group.pset>

Additional Info:
    First group defined in <group.pset> will be used
    as control group and the second group will be
    used as case group. Other groups will be ignored.

Output:
    <dataset>.asca : association analysis based on alleles
");
    abort_program();
};

string filePrefix = "";
string setFile = "";

int main(string[] args){
    flagsNoArg["-h"]  = &help;
    flags["-g"] =
        function void(int idx, string value){ setFile = value; };
    flags["def"] =
        function void(int idx, string value){
            if( filePrefix == "")
                filePrefix = value;
            else{
                writefln("Invalid argument %s", value);
                help(0,"");
                abort_program();
            }
        };
    parseFlags(args);

    if((filePrefix == "") || (setFile == "")){
        writefln("No file specified.");
        help(0,"");
    }

    GMap map = new GMap();
    map.load(filePrefix);
    auto stat = new GMapStatistics(map);

    auto mask = maskFromFile(setFile, map);
    assert(mask.group_count >= 2, "Not enough groups in file.");

    auto count = map.countAllMarkers(mask);    

    // analysis of alleles
    auto asca = new BufferedFile( filePrefix ~ ".asca", FileMode.OutNew );
    scope(exit){asca.close();}
    
    asca.writefln("MARKER\tCHISQ\tP\t%s\t%s", stat.alleles, stat.alleles);
    for (uint i = 0; i < map.markers.count; i++){
        auto controls = stat.alleleCount(count[i][0]);
        auto controls_f = countToFrequency(controls);

        auto cases = stat.alleleCount(count[i][1]);
        auto cases_f = countToFrequency(cases);

        int degreesOfFreedom;
        auto chi = caseControlChiSquare(cases, controls, degreesOfFreedom);
        auto p = pFromChiSquare(chi, degreesOfFreedom);
        
        asca.writefln("%s\t%.4f\t%.4f\t%.4s\t%.4s", map.markers.names[i], chi, p, controls_f, cases_f);
    }    
    return 1;
}