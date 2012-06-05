module gmapassocsimp;

import std.file;
import std.stream;
import std.string;

import std.stdio;
import std.math;

import gmap.gmap;
import gmap.mask;
import gmap.statistics;
import gmap.utils;

int main(string[] args){
    GMap map = new GMap();
    map.load("andmed");
    auto stat = new GMapStatistics(map);

    auto filter = maskFromFile("grupid.pset", map);
    assert(filter.group_count >= 2, "Ei ole piisavalt gruppe.");

    auto sagedused = map.countAllMarkers(filter);

    auto asca = new BufferedFile( "v2ljund.asca", FileMode.OutNew );
    scope(exit){asca.close();}

    asca.writefln("MARKER\tCHISQ\tP\t%s\t%s", stat.alleles, stat.alleles);
    for (uint i = 0; i < map.markers.count; i++){
        auto kontroll = stat.alleleCount(sagedused[i][0]);
        auto kontroll_f = countToFrequency(kontroll);

        auto juhtum = stat.alleleCount(sagedused[i][1]);
        auto juhtum_f = countToFrequency(juhtum);

        int vabadusAste;
        auto chi = caseControlChiSquare(juhtum, kontroll, vabadusAste);
        auto p = pFromChiSquare(chi, vabadusAste);

        asca.writefln("%s\t%.4f\t%.4f\t%.4s\t%.4s",
                        map.markers.names[i], chi, p, kontroll_f, juhtum_f);
    }
    return 1;
}