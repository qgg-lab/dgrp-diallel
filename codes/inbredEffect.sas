data diallel;
    infile 'diallelWolbaAdjusted.csv' dlm = ',';
    input FemaleParent $ MaleParent $ NTotal Adjusted inbred $;
run;

title 'inbred effect';
proc mixed data = diallel covtest;
    class FemaleParent MaleParent inbred;
    model Adjusted = inbred / solution;
    random FemaleParent*MaleParent(inbred);
run;

