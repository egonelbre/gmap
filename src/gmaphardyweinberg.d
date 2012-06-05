module gmaphardyweinberg;

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
"Usage: gmaphardyweinberg <dataset>

Program does a basic Hardy-Weinberg Equilibrium analysis.

Input:
    There must be files <dataset>.gmap, <dataset>.gidx,
    <dataset>.midx and <dataset>.aidx

Output:
    <dataset>.hwb : results
");
    abort_program();
};

string filePrefix = "";

int main(string[] args)
{
    flagsNoArg["-h"]  = &help;
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
    auto markers = map.countAllMarkers();

    auto hwb = new BufferedFile(filePrefix ~ ".hwb", FileMode.OutNew);
    scope(exit){ hwb.close(); }

    auto stat = new GMapStatistics(map);

    hwb.writefln("SNP\tCHI\tp\tGT\tGT(c)");
    for (uint i = 0; i < markers.length; i++){
        ulong total_alleles = 2 * sumarr(markers[i][1..$]);
        if ( total_alleles == 0){
            continue; // no valuable data at that marker
        }
            
        // calculate allele percentages
        auto p = stat.allelePercentage(markers[i]);

        auto dif_allelecount = 0;
        for(ulong j = 0; j < p.length; j++)
            if (p[j] > 0)
                dif_allelecount++;
        
        // calculate genotype expected values
        real[] expected = new real[map.genotypes.count];
        expected[] = total_alleles / 2;
        for (uint j = 1; j < map.genotypes.count; j++){
            if (!stat.homozygote[j]) // for AB and BA
                expected[j] *= 2;
            for (uint k = 0; k < p.length; k++)
                expected[j] *= pow( p[k], cast(uint)stat.gtmul[k][j] );
        }
        
        // calculate the Pearson Chi-square test
        real chi = 0;
        for (uint j = 1; j < map.genotypes.count; j++){
            // j = 1 because we ignore 0 0
            if (expected[j] == 0)
                continue; // expected[j] can only be 0 if 
                
            real diff = abs(markers[i][j] - expected[j]);
            
            // we ignore those where expected is < 5
            // because they generate Type I errors
            if (expected[j] >= 5)
                chi += diff*diff / expected[j];
        }
        chi = chi;
        
        real pvalue;
        switch ( dif_allelecount ){
            case 1,2: pvalue = pFromChiSquare(chi, 1); break;
            case   3: pvalue = pFromChiSquare(chi, 3); break;
            default : pvalue = pFromChiSquare(chi, 6); break;
        }
        
        hwb.writef("%s\t%.4f\t%.4f\t", map.markers.names[i], chi, pvalue);
        
        if (!map.markers.mm_mode){
            string genotypes = "";
            string genotypecount = "";
            for (uint j = 1; j < map.genotypes.count; j++){
                if (markers[i][j] > 0){
                    if(genotypes == ""){
                        genotypes = map.genotypes.names[j];
                        genotypecount = format("%d", markers[i][j]);
                    } else {
                        genotypes ~= "," ~ map.genotypes.names[j];
                        genotypecount ~= format(",%d", markers[i][j]);
                    }
                }
            }
            hwb.writef("[%s]\t[%s]\n", genotypes, genotypecount);
        }else{
            string fixName(string inp){ return tr(inp, "AB", map.markers.major[i] ~ map.markers.minor[i]); }
            hwb.writef("[%s,%s,%s]\t", fixName(map.genotypes.names[1]),
                                        fixName(map.genotypes.names[2]),
                                        fixName(map.genotypes.names[3]));
            hwb.writef("[%d,%d,%d]\n", markers[i][1], markers[i][2], markers[i][3]);
        }        
    }
    return 1;
}