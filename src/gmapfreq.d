module gmapfrqa;

import std.file;
import std.stream;
import std.string;

import std.stdio;
import std.math;
import std.perf;

import gmap.gmap;
import gmap.mask;
import gmap.utils;
import gmap.statistics;

import tools.flags;
import tools.bitutils;

void help(int idx, string value) {
    writef(
"Usage: gmapfreq (-c) (-g <group.pset>) <dataset>

Program counts frequencies of genotypes and alleles from
a dataset.

Input:
    There must be files <dataset>.gmap, <dataset>.gidx,
    <dataset>.midx, <dataset>.aidx

Additional Parameters:
    -c :instead of frequencies output count of genotypes
        and alleles
    -g :process for each group specified in <group.pset>

Output:
    If there were no groups specified:
        <dataset>.frqa : for alleles
        <dataset>.frqg : for genotypes
    If there was a group file specified
        <dataset>.<groupname>.frqa : for alleles
        <dataset>.<groupname>.frqg : for genotypes
");
    abort_program();
};

bool   doCount = false;
string filePrefix = "";
string groupFilename = "";

int main(string[] args){   
    flagsNoArg["-h"]  = &help;
    flags["-g"] =
        function void(int idx, string value){ groupFilename = value; };
    flagsNoArg["-c"] =
        function void(int idx, string value){ doCount = true; };
    flags["def"] =
        function void(int idx, string value){
            if( filePrefix == "")
                filePrefix = value;
            else{
                writefln(filePrefix);
                writefln("Invalid argument %s", value);
                help(0,"");
            }
        };
    parseFlags(args);

    if(filePrefix == ""){
        writefln("No file specified.");
        help(0,"");
    }

    GMap map = new GMap();
    map.load(filePrefix);

    // choose right output function
    void function(BufferedFile f, GMapStatistics stat, ulong[] count) alleleOutput, genotypeOutput;
    if (!doCount){
        alleleOutput =
            function void (BufferedFile f, GMapStatistics stat, ulong[] count){
                f.writefln(stat.allelePercentage(count));
            };
        genotypeOutput =
            function void (BufferedFile f, GMapStatistics stat, ulong[] count){
                f.writefln(stat.genotypePercentage(count));
            };
    } else {
        alleleOutput =
            function void (BufferedFile f, GMapStatistics stat, ulong[] count){
                f.writefln(stat.alleleCount(count));
            };
        genotypeOutput =
            function void (BufferedFile f, GMapStatistics stat, ulong[] count){
                f.writefln(count);
            };
    }
    
    auto stat = new GMapStatistics(map);
    if (groupFilename == ""){
        // count everything
        auto count = map.countAllMarkers();
        
        auto frqa  = new BufferedFile(filePrefix ~ ".frqa", FileMode.OutNew);
        scope(exit){frqa.close();}
        // output
        frqa.writefln("CHR\tSNP\t", stat.alleles);
        for (ulong i = 0; i < map.markers.count; i++){
            frqa.writef("%s\t%s\t", map.markers.chromosomes[i], map.markers.names[i]);
            alleleOutput(frqa, stat, count[i]);
        }

        auto frqg  = new BufferedFile(filePrefix ~ ".frqg", FileMode.OutNew);
        scope(exit){frqg.close();}

        frqg.writefln("CHR\tSNP\t", stat.genotypes);
        for (ulong i = 0; i < map.markers.count; i++){
            frqg.writef("%s\t%s\t", map.markers.chromosomes[i] , map.markers.names[i]);
            genotypeOutput(frqg, stat, count[i]);
        }
    } else {
        auto filter = maskFromFile(groupFilename, map);
        
        // count everything
        auto count = map.countAllMarkers(filter);
        
        for(int k = 0; k < filter.group_count; k++){
            auto frqa  = new BufferedFile(filePrefix ~ "." ~ filter.names[k] ~ ".frqa", FileMode.OutNew);
            scope(exit){frqa.close();}
            // output
            frqa.writefln("CHR\tSNP\t", stat.alleles);
            for (ulong i = 0; i < map.markers.count; i++){
                frqa.writef("%s\t%s\t", map.markers.chromosomes[i], map.markers.names[i]);
                alleleOutput(frqa, stat, count[i][k]);
            }
            
            auto frqg  = new BufferedFile(filePrefix ~ "." ~ filter.names[k] ~ ".frqg", FileMode.OutNew);
            scope(exit){frqg.close();}
            
            frqg.writefln("CHR\tSNP\t", stat.genotypes);
            for (ulong i = 0; i < map.markers.count; i++){
                frqg.writef("%s\t%s\t", map.markers.chromosomes[i] , map.markers.names[i]);
                genotypeOutput(frqg, stat, count[i][k]);
            }
        }
    }
    return 1;
}