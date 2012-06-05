module gmapmulassoc;

import std.file;
import std.stream;
import std.string;
import std.conv;

import std.stdio;
import std.math;
import std.random;

import std.perf;

import gmap.gmap;
import gmap.mask;
import gmap.statistics;
import gmap.utils;

import tools.flags;

void help(int idx, string value) {
    writef(
"Usage: gmapmulassoc (-c TESTS) (-b CASES) <dataset>

Program does multiple random basic case/control association
analysis to calculate baseline for randomized trial
noise.

Input:
    There must be files <dataset>.gmap, <dataset>.gidx,
    <dataset>.midx, <dataset>.aidx

Arguments:
    -c :how many case/control tests to run (default 1000)
    -b :percentage of people to be included in cases (default 20%%)

Output:
    <dataset>.asca : association analysis based on alleles
");
    abort_program();
};

string filePrefix = "";
int    testCount = 1000;
int    histogramCount = 100;
int    casePercentage = 20;

int main(string[] args){
    flagsNoArg["-h"]  = &help;
    flags["-c"] =
        function void(int idx, string value){ testCount = toInt(value); };
    flags["-h"] =
        function void(int idx, string value){ histogramCount = toInt(value); };
    flags["-b"] =
        function void(int idx, string value){ casePercentage = limits(toInt(value),0,100); };
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

    if(filePrefix == ""){
        writefln("No file specified.");
        help(0,"");
    }

    GMap map = new GMap();
    map.load(filePrefix);
    auto stat = new GMapStatistics(map);
    ulong[] masklist = new ulong[map.people.count];
    for(uint k = 0; k < masklist.length; k++)
        masklist[k] = k;
    uint pivot = masklist.length * casePercentage / 100;
    
    ulong[] maskCase;
    ulong[] maskControl;
    real[] minPValues = new real[testCount];
//     ulong[] pValueCount = new ulong[histogramCount];
//     pValueCount[] = 0;

//     real base = pow(pow(10.0,-10.0),1.0/histogramCount);
//     real pValueFromHash(uint value){
        //return 0.5*cos(value*PI/histogramCount) + 0.5;
        //return cast(real) value/cast(real) histogramCount;
//         return pow(base, value);
//     };
    
//     uint hashPValueTo(real value){
        //uint val = cast(uint) (histogramCount/PI * acos(2*(value-0.5)));
        //uint val = cast(uint) (value * histogramCount);
//         uint val = cast(uint) (log(value) / log(base));
//         if (val < 0) val = 0;
//         else if (val >= histogramCount)
//             val = histogramCount - 1;
//         return val;
//     };

    auto log = new BufferedFile("pvalues.log", FileMode.Append);
    scope(exit){log.close();}
    
    for( uint testNr = 0; testNr < testCount; testNr++){        
        randomizeArray(masklist);
        maskCase = masklist[0..pivot].sort;
        maskControl = masklist[pivot..$].sort;
        
        auto mControl = new GMapMaskSingle(maskControl, "Control");
        auto mCase = new GMapMaskSingle(maskCase, "Case");
        auto countControl = map.countAllMarkers(mControl);
        auto countCase = map.countAllMarkers(mCase);
        
        real minP = 1;
        ulong idx = 0;

        for (uint i = 0; i < map.markers.count; i++){
            auto controls = stat.alleleCount(countControl[i]);
            auto cases = stat.alleleCount(countCase[i]);
            int degreesOfFreedom;
            auto chi = caseControlChiSquare(cases, controls, degreesOfFreedom);
            auto p = pFromChiSquare(chi, degreesOfFreedom);
//             pValueCount[ hashPValueTo(p) ]++;
            
            if ( p < minP ){
                minP = p;
                idx = i;
            }
        }
        minPValues[testNr] = minP;
        writefln("#%d\t%.16f\t%s\t%s",testNr, minP, countControl[idx], countCase[idx]);
        log.writefln("%.16f\t%s\t%s", minP, countControl[idx], countCase[idx]);
        log.flush();
    }

//     auto hist = new BufferedFile("pvalues.csv", FileMode.OutNew);
//     scope(exit){hist.close();}
// 
//     hist.writef("%.6f", pValueFromHash(0));
//     for( uint i = 1; i < histogramCount; i++)
//         hist.writef(",%.10f", pValueFromHash(i));
//     hist.writefln();
//     hist.writef("%d", pValueCount[0]);
//     for( uint i = 1; i < histogramCount; i++)
//         hist.writef(",%d", pValueCount[i]);
//     hist.writefln();
    
    return 1;
}