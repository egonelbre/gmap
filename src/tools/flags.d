module tools.flags;

private import std.stdio;
private import std.c.stdlib;

void abort_program(){
    std.c.stdlib.exit(std.c.stdlib.EXIT_FAILURE);
}

alias void function (int idx, string value) flagfunc;
flagfunc[string] flags, flagsNoArg;

static this() {
    flags["err"] = function void(int idx, string value){
        abort_program();
    };
}

void parseFlags(string[] args){
    for (uint i = 1; i < args.length; i++){
        string arg = args[i];
        if(arg[0] == '-'){
            assert(arg.length >= 2, "Invalid flag " ~ arg);
            string flag = arg[0..2];
            if(flag in flags){
                string value = arg[2..$];
                if((value.length == 0) && (i+1 < args.length)){
                    value = args[i+1];
                    flags[flag](i+1,value);
                    i++;
                }
                else
                    flags[flag](i,value);                
            } else if (flag in flagsNoArg){
                flagsNoArg[flag](i,arg);
            } else {
                writefln("Invalid flag " ~ flag);
                flags["err"](i, flag);
                continue;
            }
        }
        else
            flags["def"](i,arg);
    }
}

