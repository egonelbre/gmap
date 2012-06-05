module gmaprandpheno;

import std.file;
import std.stream;
import std.string;

import std.stdio;
import std.math;
import std.perf;
import std.random;

import gmap.gmap;
import gmap.mask;
import gmap.utils;
import gmap.statistics;

import tools.flags;

void help(int idx, string value) {
    writef(
"Usage: gmaprandpheno -i <dataset> <output>

Program generates a randomized phenotype sample for doing
case/control study.

Input:
    There must be files <dataset>.gmap, <dataset>.gidx,
    <dataset>.midx, <dataset>.aidx

Output:
    <output>: columns are Family ID, Individual Id, Phenotype
    <output>.pset: pset fileformat
");
    abort_program();
};

string filePrefix = "";
string output = "";

int main(string[] args){
    flagsNoArg["-h"]  = &help;
    flags["-i"] =
        function void(int idx, string value){ filePrefix = value; };
    flags["def"] =
        function void(int idx, string value){
            if( output == "")
                output = value;
            else{
                writefln("Invalid argument %s", value);
                help(0,"");
            }
        };
    parseFlags(args);

    if(filePrefix == "" || output == ""){
        writefln("One of the files not specified.");
        help(0,"");
    }

    GMap map = new GMap();
    map.load(filePrefix);

    auto pheno = new BufferedFile(output, FileMode.OutNew);
    auto phenoset = new BufferedFile(output ~ ".pset", FileMode.OutNew);
    scope(exit){ pheno.close(); phenoset.close(); }

    ulong case_count = 0;
    ulong[] cases = new ulong[map.people.count];
    ulong control_count = 0;
    ulong[] controls = new ulong[map.people.count];

    string[] family_id;
    string[] individual_id;

    for (ulong i = 0 ; i < map.people.count; i++) {
        if (rand() % 100 == 0){
            pheno.writefln("%s\t%s\t%s", map.people.family_id[i], map.people.individual_id[i], "0" );
        }
        else if (rand() % 100 < 50){
            pheno.writefln("%s\t%s\t%s", map.people.family_id[i], map.people.individual_id[i], "2" );
            cases[case_count] = i;
            case_count++;
        } else {
            pheno.writefln("%s\t%s\t%s", map.people.family_id[i], map.people.individual_id[i], "1" );
            controls[control_count] = i;
            control_count++;
        }
    }

    phenoset.writef("control\n  ");
    for (ulong i = 0; i < control_count; i++){
        if (i % 20 == 19)
            phenoset.writef("\n  ");
        phenoset.writef("%d\t", controls[i]);
    }
    phenoset.writefln("\nEND");
    phenoset.writefln("");
    phenoset.writef("case\n  ");
    for (ulong i = 0; i < case_count; i++){
        if (i % 20 == 19)
            phenoset.writef("\n  ");
        phenoset.writef("%d\t", cases[i]);
    }
    phenoset.writefln("\nEND");

    return 1;
}