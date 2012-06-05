module gmap.mask;

import std.stdio;
import std.math;

import gmap.gmap;

enum GMapMaskTag { ADDRESS, ID, LENGTH, MASK }

/**
 * Single GMap Mask
 * Format : #addr, #len, #mask, #mask, ... 
 */
class GMapMaskSingle {
    intg[] data;
    ulong address_count = 0;
    ulong mask_count = 0;
    string name = "";

    private this(){};
    
    /***********************************
     * Constructor
     * Params:
     *    people = sorted list of people to create GMap Mask for
     */
    this ( ulong[] people, string name = "" ){
        this.name = name;
        ulong pos = people[0] / intgbits;
        address_count = 1; // count of length ids
        mask_count = 1; // count of masks
        
        ulong next_pos;
        size_t i;
        
        // count how many int masks there will be
        for (i = 1; i < people.length; i++){
            next_pos = people[i] / intgbits;
            if( pos != next_pos ){
                if( pos + 1 != next_pos )
                    address_count++;
                pos = people[i] / intgbits;
                mask_count++;
            }
        }
        // set data size
        data.length = address_count * 2 + mask_count;
        
        // populate data        
        ulong mask, data_pos, counter_pos;
        pos = people[0] / intgbits;
        data[0] = pos;
        data[1] = 0;
        counter_pos = 1;
        data_pos = 2;
        mask = 0;
        for (i = 0; i < people.length; i++){
            next_pos = people[i] / intgbits;
            if ( pos != next_pos ) {
                // add mask
                data[data_pos] = mask;
                data[counter_pos]++;
                data_pos++;
                mask = 0;
                
                // if new address should be added
                if(pos + 1 < next_pos){
                    data[data_pos] = next_pos; data_pos++;
                    data[data_pos] = 0; 
                    counter_pos = data_pos;
                    data_pos++;
                }
                pos = next_pos;
            }
            mask |= (1L << (people[i] % intgbits));
        }
        if (mask != 0 ) {
            data[data_pos] = mask;
            data[counter_pos]++;
        }
    }
}

/**
 * Multiple GMap Masks
 * Format : #addr, #id, #len, #mask, #mask, ... 
 */
class GMapMaskMultiple {
    intg[] data;
    ulong address_count = 0;
    ulong mask_count = 0;
    ulong group_count = 0;
    string[] names;

    private this(){};
    
    this (GMapMaskSingle group){
        names = new string[1];
        names[0] = group.name;
        address_count = group.address_count;
        mask_count = group.mask_count;
        group_count = 1;
        
        data.length = address_count * 3 + mask_count;
        size_t i, j; 
        ulong count;
        
        GMapMaskTag state = GMapMaskTag.ADDRESS;
        
        i,j = 0;
        while ( i < group.data.length ){
            switch ( state ){
                case GMapMaskTag.ADDRESS:
                    data[j] = group.data[i]; j++; 
                    data[j] = 0; j++; //index
                    state = GMapMaskTag.LENGTH;
                    break;
                case GMapMaskTag.LENGTH:
                    count = group.data[i];
                    data[j] = count; j++;
                    state = GMapMaskTag.MASK;
                    break;
                case GMapMaskTag.MASK:
                    data[j] = group.data[i]; j++; 
                    count--; 
                    if ( count <= 0 ) 
                        state = GMapMaskTag.ADDRESS;
                    break;
            }
            i++;
        }
    }
    
    this (GMapMaskMultiple group){
        data.length = group.data.length;
        data[] = group.data[];
        names = new string[group.names.length];
        names[] = group.names[];
        address_count = group.address_count;
        mask_count = group.mask_count;
        group_count = group.group_count;
    }
    
    this (ulong[] people, string name = ""){
        GMapMaskSingle group = new GMapMaskSingle(people, name);
        this(group);
    }
   
    GMapMaskMultiple merge(GMapMaskSingle grp){
        return merge( new GMapMaskMultiple(grp) );
    }

    GMapMaskMultiple merge(GMapMaskMultiple grpb){
        ulong addra, addrb;
        ulong addr, len;
        GMapMaskMultiple grpa = this;
        
        GMapMaskMultiple res = new GMapMaskMultiple();
        
        res.address_count = grpa.address_count + grpb.address_count;
        res.mask_count = grpa.mask_count + grpb.mask_count;
        res.group_count = grpa.group_count + grpb.group_count;
        res.data.length = grpa.data.length + grpb.data.length;

        res.names.length = grpa.names.length + grpb.names.length;
        res.names[0..grpa.names.length] = grpa.names[];
        res.names[grpa.names.length..length] = grpb.names[];
        
        addra = 0; addrb = 0;
        addr = 0;
        while( addr < res.data.length ){
            void appendGrpA(){
                res.data[addr] = grpa.data[addra]; addr++; addra++; // address
                res.data[addr] = grpa.data[addra]; addr++; addra++; // id
                len = grpa.data[addra];
                res.data[addr] = len;
                addr++; addra++; // length
                res.data[addr..addr+len] = grpa.data[addra..addra+len]; // data
                addr = addr + len; addra = addra + len;
            }
            
            void appendGrpB(){
                res.data[addr] = grpb.data[addrb]; addr++; addrb++; // address
                res.data[addr] = grpa.group_count + grpb.data[addrb]; addr++; addrb++; // id
                len = grpb.data[addrb];
                res.data[addr] = len;
                addr++; addrb++; // length
                res.data[addr..addr+len] = grpb.data[addrb..addrb+len]; // data
                addr = addr + len; addrb = addrb + len;
            }
            
            if ((addra < grpa.data.length) && (addrb < grpb.data.length)){
                if (grpa.data[addra] <= grpb.data[addrb])
                    appendGrpA();
                else 
                    appendGrpB();
            }
            else {
                if (addra < grpa.data.length)
                    appendGrpA();
                else if (addrb < grpb.data.length)
                    appendGrpB();
                else
                    break;
            }            
        }
        return res;
    }
}