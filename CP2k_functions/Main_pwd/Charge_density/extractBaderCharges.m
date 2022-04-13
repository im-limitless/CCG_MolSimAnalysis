function [Q, Qnet, StepNum] = extractBaderCharges(fldrname, flname, Atoms, AtomList)

fid = fopen([fldrname flname]);
    
    % % find the step number from the Bader file append
    StepIndx = find(flname == '_');
    StepIndx = [StepIndx find(flname == '.')];
    StepNum = str2num(flname(StepIndx(1)+1:StepIndx(2)-1));
    
    eofstat = false;
    
    % % Read the Bader charge ACF file
    HeadLine = fgetl(fid); eofstat = feof(fid); EndLine = fgetl(fid);
    i = 0;
    Bader =[];
    
    while ~eofstat
        i=i+1;
        textLine = fgetl(fid); eofstat = feof(fid);
        if strcmp(textLine, EndLine)
            break
        else
            Bader = [Bader; str2num(textLine)];
        end
    end
    
    fclose(fid);
    
    % % Extract the raw Bader charge, Q
    Q = Bader(:,5);
    
    % % Loop over the atom types/number of each type & convert raw to net charges
    TotalCharge = zeros(length(Atoms),1);
    NetCharge = zeros(length(Atoms),1);
    StndElec = zeros(length(Atoms),1);
    StockAtoms = {'Cl' 'H' 'O' 'Al'};
    Elecrons = [7; 1; 6; 3];
    nAtoms=zeros(length(AtomList),1);
    for i = 1:size(AtomList,1)
        for j = 1:size(Atoms,1)
            if  contains(AtomList(i,:), Atoms(j,:))
                TotalCharge(j,1) = Bader(j,5);
                Indx(j,i) = 1;
                for k = 1:length(StockAtoms)
                    if contains(AtomList(i,:), StockAtoms{k})
                        nAtoms(i,1) = nAtoms(i,1)+1;
                        NetCharge(j,1) = Elecrons(k,1);
                    end
                end
            end
        end
    end
    
    AlIndxs = find(ismember(AtomList, 'Al'));
    AlIndxs = AlIndxs(AlIndxs < length(AtomList)+1);
    
    NetCharge = (NetCharge-TotalCharge).*Indx;
    Qnet = sum(NetCharge,2);
      
    return