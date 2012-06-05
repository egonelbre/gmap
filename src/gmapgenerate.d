module gmapgenerate;

import std.stdio;
import std.file;
import std.math;
import std.stream;
import std.string;
import std.conv;
import std.random;

import tools.flags;

int marker_count = 1000;
int people_count = 100;

string filePrefix = "";
string genotypeFile = "";

void help(int idx, string value) {
    writef(
"Usage: gmapgenerate (-s MARKERSxPEOPLE) -g <genotypes.gt> <ouput>

Program generates random .ped and .map file. They can only
be used for speed testing as their data contents is not of
real value. Genotypes used will be picked from <genotypes.gt>.

Input:
    There must be file <genotypes.gt>

Arguments:
    -s :dimension of generated file (default 1000x100)
        expected file size will be:
            3.82MB * MARKERS * PEOPLE / 10^6
            
        5000x1000   ~   19MB
        80000x2000  ~  611MB
        200000x2000 ~ 1528MB

        Program doesn't create files larger than 2GB.
        To get larger files append two generated files.
        
        for testing it should be at least 1GB otherwise
        times will be very small
        
    -g :genotypes file
    -e :error rate (how many '0 0' will be generated)
        

Output:
    <dataset>.ped : genotype information
    <dataset>.map : marker information
");
    abort_program();
};

int main(string[] args)
{
    flagsNoArg["-h"] = &help;
    flags["-s"] =
        function void(int idx, string value){
            string[] dim = split(value,"x");
            assert(dim.length == 2, "dimensions badly specified");
            marker_count = toInt(dim[0]);
            people_count = toInt(dim[1]);
        };
    flags["-g"] =
        function void(int idx, string value){ genotypeFile = value; };
    flags["def"] =
        function void(int idx, string value){ filePrefix = value; };

    parseFlags(args);

    if (filePrefix == "" || genotypeFile == ""){
        writefln("No output or genotype file specified!");
        help(0,"");
        abort_program();
    }

    string pedfile = filePrefix ~ ".ped";
    string mapfile = filePrefix ~ ".map";

    // input files
    auto all = new BufferedFile(genotypeFile);
    scope(exit){all.close();}
    
    // output files
    auto ped = new BufferedFile(pedfile, FileMode.OutNew);
    auto map = new BufferedFile(mapfile, FileMode.OutNew);
    scope(exit){
        ped.close();
        map.close();
    }

    real approx_size = (people_count / 1000.0) * (marker_count / 1000.0) * 3.82;
    writefln("Approximate filesize will be %.2fMB",  approx_size);
    if (approx_size > 1800){
        writefln("File is too big!");
        abort_program();
    }
    // create to ".gt" mapping for genotype -> idx
    writefln("reading genotypes...");
    string[uint] gthash;
    bool[char] allelehash;
    int genotypecount = 0;
    while( !all.eof() ){ 
        string genotype = all.readLine();
        foreach(c; genotype)
            allelehash[c] = true;
        string[] genotypes = split(genotype, "\t");
        for(int i = 0; i < genotypes.length; i++)
            gthash[genotypecount] = genotypes[i];
        genotypecount++;
    }

    string alleles = "";
    foreach(allele; allelehash.keys)
        if( allele != '\t' && allele != ' ' && allele != '0' )
            alleles ~= allele;
    
    string[][] snp_genos= new string[][marker_count];
    for (int m = 0; m < marker_count; m++){
        snp_genos[m] = new string[4];
        uint a1 = rand() % alleles.length;
        uint a2 = rand() % alleles.length;
        if( a1 == a2 )
            a2 = (a2 + rand() % (alleles.length - 1)) % alleles.length;
        if( a1 == a2 )
            a2 = (a2 + 1)  % alleles.length;
        snp_genos[m][0] = alleles[a1] ~ " " ~ alleles[a1];
        snp_genos[m][1] = alleles[a1] ~ " " ~ alleles[a2];
        snp_genos[m][2] = alleles[a2] ~ " " ~ alleles[a1];
        snp_genos[m][3] = alleles[a2] ~ " " ~ alleles[a2];
    }
    
    writefln("generating data...");
    for (int p = 0; p < people_count; p++){
        string platename = format("A%d#Plate%d#%d",p,p,p);
        ped.writef("%s\t%s\t0\t0\t0\t0", platename, platename);
        for (int m = 0; m < marker_count; m++){
            if (rand() % 40 == 1)
                ped.writef("\t0 0");
            else {
                uint n = rand() % alleles.length;
                ped.writef("\t%s", snp_genos[m][n]);
            }
        }
        ped.writefln();
        ped.flush();
    }

    ulong n = 0;
    for(int m = 0; m < marker_count; m++){
        n += rand() % 3;
        string marker = format("1\tMark%d\t0\t%d", m, n);
        map.writefln("%s",marker);
    }
    
    return 1;
}