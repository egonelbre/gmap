module gmap.gmap;

import std.stdio;
import std.file;
import std.stream;
import std.string;
import std.math;
import std.mmfile;
import std.conv;

import tools.bitutils;

import gmap.mask;
import gmap.utils;

public const sufMarkers = ".midx";
public const sufGenotypes = ".gidx";
public const sufPeople  = ".pidx";
public const sufGMap    = ".gmap";

alias ulong intg;
const intgbits = intg.sizeof * 8;

class Markers {
    ulong count;
    string[] chromosomes;
    string[] names;
    ulong[] genetic_distance;
    ulong[] basepair_position;
    ulong[] positions;
    string[] major;
    string[] minor;
    bool mm_mode = false; // major-minor allele mode

    this (string fileprefix) {
        File midx = new File(fileprefix ~ sufMarkers, FileMode.In);
        scope(exit){ midx.close(); }
        count = toInt(midx.readLine());
        chromosomes = new string[count];
        names = new string[count];
        genetic_distance = new ulong[count];
        basepair_position = new ulong[count];
        positions = new ulong[count];
        
        bool doTest = true;
        for (ulong i = 0; i < count; i++){
            string[] line = split(midx.readLine(), "\t");
            if (doTest){
                if (line.length == 7){
                    major = new string[count];
                    minor = new string[count];
                    mm_mode = true;
                }
                doTest = false;
            }
            chromosomes[i] = line[0];
            names[i] = line[1];
            genetic_distance[i] = toInt(line[2]);
            basepair_position[i] = toInt(line[3]);
            //division corrects from .gmap file position to
            //intg[] array position
            positions[i] = toInt(line[4]) / (intg.sizeof);
            if (mm_mode){
                major[i] = line[5];
                minor[i] = line[6];
            }
        }
    }
}

class Genotypes {
    ulong count;
    string[] names;
    
    this (string fileprefix) {
        File gidx = new File(fileprefix ~ sufGenotypes, FileMode.In);
        // this shows how many different genotypes are
        gidx.readf("%d", &count);
        gidx.readLine();
        
        names = new string[count];
        names[] = "";
        while( !gidx.eof() ){
            string[] values = split(gidx.readLine(), "\t");
            int idx = atoi(values[1]);
            names[idx] = values[0];
        }
        gidx.close();
    }
}

class People {
    ulong count;
    string[] family_id;
    string[] individual_id;
    string[] paternal_id;
    string[] maternal_id;
    string[] sex;
    string[] phenotype;
     
    this (string fileprefix) {
        File pidx = new File(fileprefix ~ sufPeople, FileMode.In);
        scope(exit){pidx.close();}
        count = atoi(pidx.readLine());
        
        family_id = new string[count];
        individual_id = new string[count];
        paternal_id = new string[count];
        maternal_id = new string[count];
        sex = new string[count];
        phenotype = new string[count];
        
        for(ulong i=0; i < count; i++){
            string[] line = split(pidx.readLine(), "\t");
            family_id[i] = line[1];
            individual_id[i] = line[2];
            paternal_id[i] = line[3];
            maternal_id[i] = line[4];
            sex[i] = line[5];
            phenotype[i] = line[6];
        }
    }
}

class GMap {
    intg[] data;
    string fileprefix;
    ulong  position;
    ulong  line_length;
    
    Markers markers;
    Genotypes genotypes;
    People  people;

    MmFile datafile;
    
    this () {}
    
    void load (string fileprefix) {
        this.fileprefix = fileprefix;
        datafile = new MmFile(fileprefix ~ sufGMap);
        data = cast(intg[]) datafile[];
        
        markers = new Markers(fileprefix);
        genotypes = new Genotypes(fileprefix);
        people  = new People(fileprefix);
        
        line_length = people.count / intgbits;
        if( people.count % intgbits != 0)
            line_length++;
            
        position = 0;
    }
    
    // counts from position for all genotypes
    ulong[] countMarkerFrom(ulong start){
        ulong[] countArr = new ulong[genotypes.count];
        ulong addr = start;
        ulong end = start + genotypes.count * line_length;

        for (auto k = 0; k < genotypes.count; k++) {
            ulong count = 0;
            end = addr + line_length;
            while (addr < end){
                count += popcount(data[addr]);
                addr++;
            }
            countArr[k] = count;
        }
        return countArr;
    }
    
    ulong[] countMarkerFrom(ulong start, GMapMaskSingle mask){
        // mask format : #addr, #len, #mask, #mask, ... 
        ulong[] countArr = new ulong[genotypes.count];
        ulong addr = start;
        ulong mask_length = mask.data.length;
        
        for (auto k = 0; k < genotypes.count; k++) {
            ulong count = 0;
            ulong mi = 0; // mask index
            ulong end;
            
            start = addr; // snp+allele start
            while (mi < mask_length){
                addr = start + mask.data[mi]; mi++;
                end = addr + mask.data[mi]; mi++;
                for (; addr < end; addr++){
                    count += popcount(data[addr] & mask.data[mi]);
                    mi++;
                }
            }
            addr = start + line_length;
            countArr[k] = count;
        }
        
        return countArr;
    }

    // result data[group][genotype]
    ulong[][] countMarkerFrom(ulong start, GMapMaskMultiple mask){
        // mask format : #addr, #id, #len, #mask, #mask, ... 
        ulong[][] countArr = new ulong[][genotypes.count]; // result
        
        ulong addr = start;
        ulong mask_length = mask.data.length;
        
        for (auto k = 0; k < genotypes.count; k++) {
            ulong mi = 0; // mask index
            ulong end, idx, icount;
            ulong[] count = new ulong[mask.group_count];
            count[] = 0;
            
            start = addr; // snp+allele start
            while (mi < mask_length){
                addr = start + mask.data[mi]; mi++;
                idx = mask.data[mi]; mi++;
                end = addr + mask.data[mi]; mi++;
                icount = count[idx];
                for (; addr < end; addr++){
                    icount += popcount(data[addr] & mask.data[mi]);
                    mi++;
                }
                count[idx] = icount;
            }
            addr = start + line_length;
            countArr[k] = count;
        }
        
        return rotateArray(countArr);
    }
    
    // counts from position for all genotypes
    // basically for performance testing
    // result data[genotype]
    ulong[] countAll(){
        ulong[] countArr = new ulong[genotypes.count];
        ulong addr = 0;
        
        auto dlen = data.length;
        ulong k = 0;
        ulong end;
        while(addr < dlen){
            ulong count = countArr[k];
            end = addr + line_length;
            while (addr < end){
                count += popcount(data[addr]);
                addr++;
            }
            countArr[k] = count;
            k = (k + 1) % genotypes.count;
        }
        return countArr;
    }

    // result data[marker][genotype]
    ulong[][] countAllMarkers(){

        ulong[][] countArr = new ulong[][markers.count];
        for(ulong i = 0; i < markers.count; i++){
            countArr[i] = countMarkerFrom(markers.positions[i]);
        }
        return countArr;
    }

    // result data[marker][genotype]
    ulong[][] countAllMarkers(GMapMaskSingle mask){
        ulong[][] countArr = new ulong[][markers.count];
        for(ulong i = 0; i < markers.count; i++){
            countArr[i] = countMarkerFrom(markers.positions[i], mask);
        }
        return countArr;
    }

    // result data[marker][group][genotype]
    ulong[][][] countAllMarkers(GMapMaskMultiple mask){
        ulong[][][] countArr = new ulong[][][markers.count];
        for(ulong i = 0; i < markers.count; i++){
            countArr[i] = countMarkerFrom(markers.positions[i], mask);
        }
        return countArr;
    }
}
