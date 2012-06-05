module gmappack;

import std.file;
import std.stream;
import std.string;

import std.stdio;
import std.math;

import gmap.gmap;
import gmap.statistics;
import gmap.utils;

import tools.flags;

void help(int idx, string value) {
    writef(
"Usage: gmappack -i <input> <output>

Program keeps SNPs from <input> that have only two
alleles all others will be discarded. Also it converts
those two alleles to A and B (A is major-allel and B
is minor-allele). Real allele names will be appended
to .midx file (first being A and second B).

Input:
    There must be files <input>.gmap, <input>.gidx,
    <input>.midx, <input>.aidx and <input.pset>

Additional Info:
    <input> and <output> can't be the same.
    Any marker discarded will be written to <output>.log.

Output:
    <output>.gmap : binary genotype data
    <output>.gidx : genotypes index
    <output>.midx : marker index
    <output>.pidx : people index
");
    abort_program();
};

string prefix = "";
string prefixout = "";

int main(string[] args)
{
    flagsNoArg["-h"]  = &help;
    flags["-i"] =
        function void(int idx, string value){ prefix = value; };
    flags["def"] =
        function void(int idx, string value){
            if( prefixout == "")
                prefixout = value;
            else{
                writefln("Invalid argument %s", value);
                help(0,"");
            }
        };
    parseFlags(args);

    if((prefix == "") || (prefixout == "")){
        writefln("No file specified.");
        help(0,"");
    }

    if(prefix == prefixout){
        writefln("Input and output can't be same.");
        help(0,"");
    }

    // log file
    auto log = new BufferedFile(prefixout ~ ".log", FileMode.OutNew);
    scope(exit){ log.close(); }

    auto map = new GMap();
    map.load(prefix);
    if (map.markers.mm_mode) {
        writefln("Data already in major-minor mode.");
        abort_program();
    }
    auto stat = new GMapStatistics(map);
    
    // our output files
    auto gmap  = new BufferedFile(prefixout ~ sufGMap, FileMode.OutNew);
    auto pidxo = new BufferedFile(prefix ~ sufPeople, FileMode.In);
    auto pidx  = new BufferedFile(prefixout ~ sufPeople, FileMode.OutNew);
    auto midx  = new BufferedFile(prefixout ~ sufMarkers, FileMode.OutNew);
    auto gidx  = new BufferedFile(prefixout ~ sufGenotypes, FileMode.OutNew);
    scope(exit){
        gmap.close();
        pidxo.close();
        pidx.close();
        midx.close();
        gidx.close();
    }

    string pline = "";
    while ( !pidxo.eof()){
        pline = pidxo.readLine();
        pidx.writefln(pline);
    }

    // generate new genotype mapping
    gidx.writefln("4");
    gidx.writefln("0 0\t0");
    gidx.writefln("A A\t1");
    gidx.writefln("A B\t2");
    gidx.writefln("B A\t2");
    gidx.writefln("B B\t3");

    ulong markercount = 0;
    string[] midxlines = new string[map.markers.count];
    
    // generate new gmap
    for(ulong i = 0; i < map.markers.count; i++){
        auto count = map.countMarkerFrom(map.markers.positions[i]);
        auto p = stat.allelePercentage(count);
        auto order = alleleOrder(p);
        if( order.length != 2){
            log.writefln(map.markers.names[i], " ", order," ", count, " ", p ," ",order.length);
            continue;
        }

        midxlines[markercount] =
            format("%s\t%s\t%d\t%d\t%d\t%s\t%s",
                map.markers.chromosomes[i],
                map.markers.names[i],
                map.markers.genetic_distance[i],
                map.markers.basepair_position[i],
                gmap.position,
                stat.alleles[order[0]],
                stat.alleles[order[1]]
            );
        markercount++;

        auto line = new intg[map.line_length];
        void copy_line(ulong base, ulong genotype){
            ulong start = base + genotype*map.line_length;
            line[] = map.data[start..start + map.line_length];
            gmap.writeBlock(line.ptr, intg.sizeof*map.line_length);
        }
        int findHomozygote(ulong a){
            for (uint i = 0; i < stat.genotypes.length; i++){
                if( stat.gtmul[a][i] == 2 )
                    return i;
            }
            assert(0, "Homozygote for " ~ stat.alleles[a] ~ stat.alleles[a] ~ " not found!");
        }
        int findHeterozygote(uint a1, uint a2){
            for (uint i = 0; i < stat.genotypes.length; i++){
                if( stat.gtmul[a1][i] == 1 && stat.gtmul[a2][i] == 1 )
                    return i;
            }
            assert(0, "Heterozygote for " ~ stat.alleles[a1] ~ stat.alleles[a2] ~ " not found!");
            return -1;
        }
        copy_line(map.markers.positions[i], 0);
        copy_line(map.markers.positions[i], findHomozygote(order[0]));
        copy_line(map.markers.positions[i], findHeterozygote(order[0], order[1]));
        copy_line(map.markers.positions[i], findHomozygote(order[1]));
    }

    midx.writefln(markercount);
    for(ulong i=0; i < markercount; i++)
        midx.writefln(midxlines[i]);
    return 1;
}
