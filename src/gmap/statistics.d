module gmap.statistics;

import std.string;
import std.stdio;
import std.math;

import gmap.gmap;
import gmap.utils;

class GMapStatistics {
    string[] alleles;
    string[] genotypes;

    ulong[][] gtmul; // shows how many alleles are in a genotype
    bool[]    homozygote; // shows whether genotypes are homozygote (true) or not (false)

    this (GMap gmap){
        genotypes = new string[gmap.genotypes.count];
        genotypes[] = gmap.genotypes.names[];
        
        string gts = join(genotypes,"");
        gts = removechars(gts, " 0");
        while( gts.length > 0 ){
            alleles ~= gts[0..1];
            gts = removechars(gts, gts[0..1]);
        }
        
        gtmul = new ulong[][alleles.length];
        for(ulong i = 0; i < alleles.length; i++){
            gtmul[i] = new ulong[genotypes.length];
            for (uint j = 0; j < genotypes.length; j++)
                gtmul[i][j] = countchars(genotypes[j], alleles[i]);
        }
        
        homozygote = new bool[genotypes.length];
        HOMOZYGOTE:
        for(uint i = 0; i < genotypes.length; i++){
            for (uint j = 0; j < alleles.length; j++){
                switch(gtmul[j][i]){
                    case 1 : homozygote[i] = false; continue HOMOZYGOTE;
                    case 2 : homozygote[i] = true; continue HOMOZYGOTE;
                    default: break;
                }
            }
        }
    }

    ulong[] alleleCount(ulong[] gtcount){
        auto res = new ulong[alleles.length];        
        for (uint j = 0; j < res.length; j++)
            res[j] = summularr(gtmul[j], gtcount);
        return res;
    }

    real[] allelePercentage(ulong[] gtcount){
        auto res = new real[alleles.length];
        auto all_count = alleleCount(gtcount);
        auto total = cast(real) sumarr(all_count);
        for (uint j = 0; j < res.length; j++)
            res[j] = cast(real) all_count[j] / total;
        return res;
    }

    real[] genotypePercentage(ulong[] gtcount){
        auto res = new real[genotypes.length];
        auto total = cast(real) sumarr(gtcount);
        for (uint j = 0; j < res.length; j++)
            res[j] = cast(real) gtcount[j] / total;
        return res;
    }
}

real pearsonChiSquare(ulong[] observed, ulong[] expected){
    assert(observed.length == expected.length, "arrays not same size");
    real res = 0;
    for (uint i = 0; i < observed.length; i++){
        if (expected[i] > 0)
            res += sqr(observed[i] - expected[i]) / cast(real)expected[i];
    }
    return res;
}

real pearsonChiSquare(real[] observed, real[] expected){
    assert(observed.length == expected.length, "arrays not same size");
    real res = 0;
    for (uint i = 0; i < observed.length; i++){
        if (expected[i] > 0)
            res += sqr(observed[i] - expected[i]) / cast(real)expected[i];
    }
    return res;
}

real caseControlChiSquare(ulong[] cases, ulong[] controls, out int degreesOfFreedom){
    assert(cases.length == controls.length, "arrays not same size");
    auto count = new ulong[cases.length];
    count[] = 0;
    degreesOfFreedom = 0;
    for(uint i = 0; i < count.length; i++){
        count[i] = cases[i] + controls[i];
        if( count[i] > 0)
            degreesOfFreedom++;
    }
    degreesOfFreedom = (degreesOfFreedom - 1) * (degreesOfFreedom - 1);
    
    ulong total_cases = sumarr(cases);
    ulong total_controls = sumarr(controls);
    ulong total = total_cases + total_controls;

    real frequency_cases = cast(real) total_cases / cast(real) total;
    real frequency_controls = cast(real) total_controls / cast(real) total;

    // expected_case[allele] =  total_cases * (count[allele] / total);
    auto expected_cases = new real[count.length];
    auto expected_controls = new real[count.length];
    for(uint i = 0; i < count.length; i++){
        expected_cases[i] = count[i] * frequency_cases;
        expected_controls[i] = count[i] * frequency_controls;
    }
    
    real res = 0;
    for (uint i = 0; i < count.length; i++)
        if (expected_cases[i] >= 5)
            res += sqr(cases[i] - expected_cases[i]) / expected_cases[i];
    for (uint i = 0; i < count.length; i++)
        if (expected_controls[i] >= 5)
            res += sqr(controls[i] - expected_controls[i]) / expected_controls[i];
    
    return res;
}

real yatesChiSquare(ulong[] observed, ulong[] expected){
    assert(observed.length == expected.length, "arrays not same size");
    real res = 0;
    for (uint i = 0; i < observed.length; i++){
        if (expected[i] > 0)
            res += sqr(abs(cast(real)observed[i] - expected[i]) - 0.5) / cast(real)expected[i];
    }    
    return res;
}

real[] countToFrequency(ulong[] count){
    real[] res = new real[count.length];
    ulong sum = sumarr(count);
    for(ulong i = 0 ; i < res.length; i++)
        res[i] = cast(real)count[i] / cast(real)sum;
    return res;
}

real pFromChiSquare(real x, int n){
    if(n==1 && x > 1000) return 0;
    if(x > 1000 || n > 1000){
        real q = pFromChiSquare((x-n)*(x-n)/(2*n),1)/2;
        if (x > n)
            return q;
        return 1-q;
    }
    real p = exp(-0.5*x);
    if (n % 2 == 1)
        p *= SQRT2*sqrt(x)/sqrt(PI);
    int k = n; while( k >= 2) { p = p * x / k; k -= 2;}
    real t = p; int a = n; while(t > 0.0000000001*p) { a += 2; t *= x / a; p += t; }
    return 1 - p;
}

// returns allele order indexes
// doesn't list alleles < 0
ulong[] alleleOrder(real[] alleles){
    ulong[] order = new ulong[alleles.length];
    for(uint k = 0; k < alleles.length; k++)
        order[k] = k;
    real[] copy = alleles.dup;
    
    bool changed = true;
    ulong last = alleles.length;
    while(changed){
        changed = false;
        last--;
        for(ulong i=0; i < last; i++){
            if (copy[i] < copy[i+1]){
                real t1 = copy[i]; copy[i] = copy[i+1]; copy[i+1] = t1;
                ulong t2 = order[i]; order[i] = order[i+1]; order[i+1] = t2;
                changed = true;
            }
        }
    }
    //trim 0 values
    last = order.length;
    while( last > 0 && copy[last-1] == 0 )
        last--;
    order.length = last;
    return order;
}