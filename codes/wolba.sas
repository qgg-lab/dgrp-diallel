data diallel;
    infile 'sas/wolbaData.csv' dlm = ',';
    input tag cross $ FemaleParent $ MaleParent $ NFemale NMale NTotal Rep $ FemaleWolba $ MaleWolba $ Wolba $;
run;

title 'Wolbachia: NTotal';
proc mixed data = diallel covtest;
    class FemaleParent MaleParent Rep FemaleWolba MaleWolba Wolba;
    model NTotal = FemaleWolba MaleWolba FemaleWolba*MaleWolba / solution;
    random FemaleParent(FemaleWolba) MaleParent(MaleWolba) FemaleParent*MaleParent(FemaleWolba*MaleWolba);
run;

title 'Wolbachia: NTotal';
proc mixed data = diallel covtest;
    class FemaleParent MaleParent Rep FemaleWolba MaleWolba Wolba;
    model NTotal = Wolba / solution;
    random FemaleParent*MaleParent(Wolba);
run;
