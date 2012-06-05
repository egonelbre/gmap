module tools.bitutils;

byte[1 << 16] bitcount;

private {
    ulong m1  = 0x5555555555555555L; //binary: 0101...
    ulong m2  = 0x3333333333333333L; //binary: 00110011..
    ulong m4  = 0x0f0f0f0f0f0f0f0fL; //binary:  4 zeros,  4 ones ...
    ulong m8  = 0x00ff00ff00ff00ffL; //binary:  8 zeros,  8 ones ...
    ulong m16 = 0x0000ffff0000ffffL; //binary: 16 zeros, 16 ones ...
    ulong m32 = 0x00000000ffffffffL; //binary: 32 zeros, 32 ones
    ulong hff = 0xffffffffffffffffL; //binary: all ones
    ulong h01 = 0x0101010101010101L; //the sum of 256 to the power of 0,1,2,3...
    ulong hff16=0xffff;
    ulong hff32=0xffffffff;
}

version(D_InlineAsm_X86){
    static uint popcount( ushort x){
        asm {
            naked                  ; // we don't want a stack frame
            and EAX, 0xFFFF        ; // just checkin
            mov AL, bitcount[EAX]  ; // get value from table
            xor AH, AH             ; // clean high byte of ax
            ret                    ; // done
        }
    }
    
    static uint popcount( uint x ){
        asm {
            naked                   ; // we don't want a stack frame
            mov ECX, EAX            ; // make a copy to EDX
            and EAX, 0xFFFF         ; // keep the lower ushort
            mov AL,  bitcount[EAX]  ; // move the lower ushort bitcount to al
            xor AH, AH;             ; // remove the upper byte of EAX
            shr ECX, 16             ; // remove the lower ushort from ebx
            add AL,  bitcount[ECX]  ; // add the upper ushort bitcount
            ret                     ; // done
        }
    }
    
    static uint popcount_b( ulong x){
        asm {
            naked;
            xor EAX, EAX            ; // clear result
            mov ECX, [ESP + 4]      ; // move ulong lower int to ecx
            mov EDX, ECX            ; // make a copy to edx
            and EDX, 0xFFFF         ; // keep the lower ushort
            add AL,  bitcount[EDX]  ; // move the lower ushort bitcount to al
            shr ECX, 16             ; // remove the lower ushort from ebx
            add AL,  bitcount[ECX]  ; // add the upper ushort bitcount
            mov ECX, [ESP + 8]      ; // move ulong higher int to ecx
            mov EDX, ECX            ; // make a copy to edx
            and EDX, 0xFFFF         ; // keep the lower ushort
            add AL,  bitcount[EDX]  ; // move the lower ushort bitcount to al
            shr ECX, 16             ; // remove the lower ushort from ebx
            add AL,  bitcount[ECX]  ; // add the upper ushort bitcount
            ret 8                   ; // clean our stack
        }
    }

    static uint popcount( ulong x){
        asm {
            naked;
            xor EAX, EAX            ; // clear result
            xor ECX, ECX            ; // clear ECX
            mov CX, [ESP + 4]       ; // get ushort
            add AL,  bitcount[ECX]  ; // add bitcount to result
            mov CX, [ESP + 6]       ; // get ushort
            add AL,  bitcount[ECX]  ; // add bitcount to result
            mov CX, [ESP + 8]       ; // get ushort
            add AL,  bitcount[ECX]  ; // add bitcount to result
            mov CX, [ESP + 10]      ; // get ushort
            add AL,  bitcount[ECX]  ; // add bitcount to result
            ret 8                   ; // clean our stack
        }
    }    

    static uint seander_popcount( ulong x){
        asm{
            naked;
            xor EAX, EAX            ; // clear result
            mov ECX, [ESP + 4]      ; // move ulong lower int to ECX
            //  v = v - ((v >> 1) & 0x55555555);
            mov EDX, ECX            ; // copy
            shr EDX, 1              ; // v >>1
            and EDX, 0x55555555     ; // ((v >>1) & 0x555555555)
            sub ECX, EDX            ; // v - ((v >>1) & 0x555555555)
            //  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
            mov EDX, ECX            ;
            shr EDX, 2              ; // v >> 2
            and EDX, 0x33333333     ; // (v >> 2) & 0x33333333
            and ECX, 0x33333333     ; // v & 0x33333333
            add ECX, EDX            ; // (v & 0x33333333) + ((v >> 2) & 0x33333333)
            //  c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
            mov EDX, ECX            ;
            shr EDX, 4              ; // v >> 4
            add ECX, EDX            ; // v + (v >> 4)
            and ECX, 0xF0F0F0F      ; // v + (v >> 4) & 0xF0F0F0F
            imul ECX, ECX, 0x01010101; // (v + (v >> 4) & 0xF0F0F0F) * 0x1010101
            shr ECX, 24             ; // ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24
            add EAX, ECX            ;
            // repeat with upper uint
            mov ECX, [ESP + 8]      ; // move ulong upper int to ECX
            //  v = v - ((v >> 1) & 0x55555555);
            mov EDX, ECX            ; // copy
            shr EDX, 1              ; // v >>1
            and EDX, 0x55555555     ; // ((v >>1) & 0x555555555)
            sub ECX, EDX            ; // v - ((v >>1) & 0x555555555)
            //  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
            mov EDX, ECX            ;
            shr EDX, 2              ; // v >> 2
            and EDX, 0x33333333     ; // (v >> 2) & 0x33333333
            and ECX, 0x33333333     ; // v & 0x33333333
            add ECX, EDX            ; // (v & 0x33333333) + ((v >> 2) & 0x33333333)
            //  c = ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
            mov EDX, ECX            ;
            shr EDX, 4              ; // v >> 4
            add ECX, EDX            ; // v + (v >> 4)
            and ECX, 0x0F0F0F0F     ; // v + (v >> 4) & 0xF0F0F0F
            imul ECX, ECX, 0x01010101; // (v + (v >> 4) & 0xF0F0F0F) * 0x1010101
            shr ECX, 24             ; // ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24
            add EAX, ECX            ;
            ret 8                   ; // clean our stack
        }
    }
}else{
    static uint popcount( ushort x){
        return bitcount[x];
    }
    
    static uint popcount( uint x ){
        return bitcount[x & 0xFFFF] +
            bitcount[(x >> 16) & 0xFFFF];
    }

    static uint popcount( ulong x ){
    return bitcount[x & 0xFFFF] +
           bitcount[(x >> 16) & 0xFFFF] +
           bitcount[(x >> 32) & 0xFFFF] +
           bitcount[(x >> 48) & 0xFFFF];
    }
}

static byte popcount_parallel( ulong x ){
    x = (x & m1 ) + ((x >>  1) & m1 );
    x = (x & m2 ) + ((x >>  2) & m2 );
    x = (x & m4 ) + ((x >>  4) & m4 );
    x = (x & m8 ) + ((x >>  8) & m8 );
    x = (x & m16) + ((x >> 16) & m16);
    x = (x & m32) + ((x >> 32) & m32);
    return x;
}

static this() {
    for(uint x = 0; x < bitcount.length; x++)
        bitcount[x] = popcount_parallel(x);
}
