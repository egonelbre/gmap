module gmap.utils;

import std.file;
import std.stream;
import std.string;
import std.conv;

import std.stdio;
import std.math;
import std.random;

import gmap.gmap;
import gmap.mask;

// calculates sum of an array
ulong sumarr(ulong[] arr){
    ulong total = 0;
    for(uint i = 0; i < arr.length; i++)
        total += arr[i];
    return total;
}

// multiplicates two arrays and then finds the sum
ulong summularr(ulong[] arr, ulong[] arr2){
    assert(arr.length == arr2.length, "arrays not same size");
    ulong res = 0;
    for (uint i = 0; i < arr.length; i++)
        res += arr[i]*arr2[i];
    return res;
}

ulong[][] rotateArray(ulong[][] arr){
    auto res = new ulong[][arr[0].length];
    for(ulong x = 0; x < res.length; x ++)
        res[x] = new ulong[arr.length];
    for(ulong x = 0; x < res.length; x ++ ){
        for(ulong y = 0; y < arr.length; y ++ ){
            res[x][y] = arr[y][x];
        }
    }
    return res;
}

void randomizeArray(ulong[] arr){
    for (uint i = 0 ; i < arr.length; i++){
        uint idx = i + rand() % (arr.length - i);
        ulong temp = arr[i];
        arr[i] = arr[idx];
        arr[idx] = temp;
    }
}

real sqr(real a){ return a * a; }

long limits(long value, long lower, long upper){
    if (value > upper)
        return upper;
    else if (value < lower)
        return lower;
    else
        return value;
}

GMapMaskMultiple maskFromFile(string filename, GMap gmap){
    GMapMaskMultiple res = null;
    string data = cast(char[]) read(filename);
    string[] elems = split(data);
    ulong[] people = new ulong[elems.length];
    ulong k = 0; // k index in people, n index in gmap.people
    string name;
    ulong state = 0; // 0 - next header, 1 - person
    for(ulong i = 0; i < elems.length; i++) {
        switch (state) {
            case 0 :
                name = elems[i];
                k = 0;
                state = 1;
                break;
            case 1 :
                if (elems[i] == "END"){
                    GMapMaskSingle grp = new GMapMaskSingle(people[0..k], name);
                    if (res !is null)
                        res = res.merge(grp);
                    else
                        res = new GMapMaskMultiple(grp);
                    k = 0;
                    state = 0;
                } else {
                    people[k] = toInt(elems[i]);
                    k++;
                }
                break;
            default : break;
        }
    }
    return res;
}