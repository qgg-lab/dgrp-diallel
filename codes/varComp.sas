data dummy;
    do row = 1 to 2500;
        output;
    end;

data fam1;
    do dam1 = 1 to 50;
        do sire1 = 1 to 50;
            output;
        end;
    end;

data fam1;
    merge dummy fam1;

data fam1;
    set fam1;
    do dam2 = 1 to 50;
        do sire2 = 1 to 50;
            output;
        end;
    end;

data dummy;
    set dummy;
    do col = 1 to 2500;
        output;
    end;

data ii;
    merge fam1 dummy;
    if row < col then delete;
    if (sire1 ^= sire2 and dam1 ^= dam2 and sire1 ^= dam2 and dam1 ^= sire2)
        then delete;   
    if sire1 = sire2 and dam1 = dam2 then rel = 'fulls';
    else if sire1 = dam2 and dam1 = sire2 then rel = 'rfs';
    else if sire1 = sire2 and dam1 ^= dam2 then rel = 'phs';
    else if sire1 ^= sire2 and dam1 = dam2 then rel = 'mhs';
    else rel = 'rhs';
    do parm = 1 to 5;
        output;
    end;

data ii;
    set ii;
    if parm = 1 then do;
        if rel = 'fulls' or rel = 'rfs' then value = 2;
        else value = 1;
        end;
    if parm = 2 then do;
        if rel = 'fulls' or rel = 'rfs' then value = 1;
        else delete;
        end;
    if parm = 3 then do;
        if rel = 'fulls' or rel = 'mhs' then value = 1;
        else delete;
        end;
    if parm = 4 then do;
        if rel = 'fulls' or rel = 'phs' then value = 1;
        else delete;
        end;
    if parm = 5 then do;
        if rel = 'fulls' then value = 1;
        else delete;
        end;
    keep parm row col value;

run;

data diallel;
    infile 'diallelWolbaAdjusted.csv' dlm = ',';
    input FemaleParent $ MaleParent $ NTotal Adjusted;
run;

title 'cross';
proc mixed data = diallel covtest;
    class FemaleParent MaleParent;
    model Adjusted = / solution;
    random FemaleParent*MaleParent;
run;

title 'Fry approach: NTotal';
proc mixed data = diallel covtest;
    class FemaleParent MaleParent;
    model Adjusted = / solution;
    random FemaleParent*MaleParent / type = lin(5) ldata = ii;
    parms / lowerb = 0,0,0,0,0,0;
run;

/* the five variance components are n nn f m fm*/
