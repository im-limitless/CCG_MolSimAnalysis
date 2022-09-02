function CreateCP2KInputfile_CP_Like(BaseIn, Base, OutFldr, Wall, CoordOut, VelOut, ABCOut, FixOut, Step, EXOR, Gam, NGam, Dir)

% write submission script


OutLoc = [Base '\' Dir];

    if exist([OutLoc],'dir')
        warning('Directory already exists! Continuing with next file.');
        return
    else
        mkdir([OutLoc]);
    end

fidout = fopen([OutLoc '\'  Dir '.inp'],'w');

fprintf(fidout,['# CP2K input file created by Matt Darby, Imperial College London for system ' OutFldr newline]);
fprintf(fidout,[' &GLOBAL' newline]);
fprintf(fidout,['   PRINT_LEVEL  LOW' newline]);
fprintf(fidout,['   PROJECT_NAME ' OutFldr newline]);
fprintf(fidout,['   RUN_TYPE  MD' newline]);
fprintf(fidout,['   WALLTIME ' num2str(Wall*60*60-300, '%i') newline]);
fprintf(fidout,[' &END GLOBAL' newline]);
fprintf(fidout,[' &MOTION' newline]);
fprintf(fidout,['   &MD' newline]);
fprintf(fidout,['     ENSEMBLE  LANGEVIN' newline]);
fprintf(fidout,['	  &LANGEVIN' newline]);
fprintf(fidout,['	   GAMMA     ' num2str(Gam) newline]);
fprintf(fidout,['	   NOISY_GAMMA     ' num2str(NGam) newline]);
fprintf(fidout,['     &END LANGEVIN' newline]);
fprintf(fidout,['     STEPS  1000' newline]);
fprintf(fidout,['     TIMESTEP     0.5' newline]);
%     fprintf(fidout,['     TEMPERATURE     3.4000000000000000E+02' newline]);
%     fprintf(fidout,['     TEMP_TOL     70.0000000000000000E+00' newline]);
%     fprintf(fidout,['	 &THERMOSTAT' newline]);
%     fprintf(fidout,['       &CSVR' newline]);
%     fprintf(fidout,['         TIMECON     1000' newline]);
%     fprintf(fidout,['       &END CSVR' newline]);
%     fprintf(fidout,['     &END THERMOSTAT' newline]);
fprintf(fidout,['    &END MD' newline]);
fprintf(fidout,['   &CONSTRAINT' newline]);
fprintf(fidout,['     &FIXED_ATOMS' newline]);
fprintf(fidout,[FixOut newline]);
fprintf(fidout,['     &END FIXED_ATOMS' newline]);
fprintf(fidout,['   &END CONSTRAINT' newline]);
fprintf(fidout,['   &PRINT' newline]);
fprintf(fidout,['     &VELOCITIES  ON' newline]);
fprintf(fidout,['     &END VELOCITIES' newline]);
fprintf(fidout,['     &RESTART  SILENT' newline]);
fprintf(fidout,['       ADD_LAST  NUMERIC' newline]);
fprintf(fidout,['       &EACH' newline]);
fprintf(fidout,['         MD  1' newline]);
fprintf(fidout,['       &END EACH' newline]);
fprintf(fidout,['     &END RESTART' newline]);
fprintf(fidout,['   &END PRINT' newline]);
fprintf(fidout,[' &END MOTION' newline]);
fprintf(fidout,[' &FORCE_EVAL' newline]);
fprintf(fidout,['   METHOD  QS' newline]);
%     fprintf(fidout,['   &EXTERNAL_POTENTIAL' newline]);
%     fprintf(fidout,['    ATOMS_LIST ' num2str(IonicSpecs') newline]);
%     fprintf(fidout,['    FUNCTION (1.0E-6)*((Z-' num2str(Vector(3,3)*1.88973/2) ')^4)' newline]);
%     fprintf(fidout,['   &END EXTERNAL_POTENTIAL' newline]);
fprintf(fidout,['   &DFT' newline]);
fprintf(fidout,['     BASIS_SET_FILE_NAME ./GTH_BASIS_SETS' newline]);
fprintf(fidout,['     POTENTIAL_FILE_NAME ./GTH_POTENTIALS' newline]);
fprintf(fidout,['     WFN_RESTART_FILE_NAME ../../' OutFldr '-RESTART.wfn' newline]);
fprintf(fidout,['     &SCF' newline]);
%     fprintf(fidout,['       MAX_SCF  800' newline]);
fprintf(fidout,['       MAX_SCF_HISTORY  3' newline]);
fprintf(fidout,['       EPS_SCF     1.0E-06' newline]);
fprintf(fidout,['       EPS_SCF_HISTORY     1.0E-05' newline]);
fprintf(fidout,['       EPS_DIIS     1.0E-01' newline]);
fprintf(fidout,['       SCF_GUESS  HISTORY_RESTART' newline]);
fprintf(fidout,['       &OT  T' newline]);
fprintf(fidout,['         MINIMIZER  DIIS' newline]);
fprintf(fidout,['         SAFE_DIIS  F' newline]);
fprintf(fidout,['         N_HISTORY_VEC  7' newline]);
fprintf(fidout,['         STEPSIZE     ' num2str(Step) newline]);
fprintf(fidout,['         PRECONDITIONER  FULL_SINGLE_INVERSE' newline]);
fprintf(fidout,['         ENERGY_GAP     1.0000000000000000E-03' newline]);
fprintf(fidout,['       &END OT' newline]);
fprintf(fidout,['       &MIXING  T' newline]);
fprintf(fidout,['         METHOD  DIRECT_P_MIXING' newline]);
fprintf(fidout,['         ALPHA     3.000E-01' newline]);
fprintf(fidout,['       &END MIXING' newline]);
fprintf(fidout,['       &PRINT' newline]);
fprintf(fidout,['         &RESTART_HISTORY  SILENT' newline]);
fprintf(fidout,['           BACKUP_COPIES  ' num2str(EXOR+2) newline]);
fprintf(fidout,['           &EACH' newline]);
fprintf(fidout,['           &END EACH' newline]);
fprintf(fidout,['         &END RESTART_HISTORY' newline]);
fprintf(fidout,['       &END PRINT' newline]);
fprintf(fidout,['     &END SCF' newline]);
fprintf(fidout,['     &QS' newline]);
fprintf(fidout,['       EPS_DEFAULT     1.0E-10' newline]);
fprintf(fidout,['       MAP_CONSISTENT  T' newline]);
fprintf(fidout,['       EXTRAPOLATION  ASPC' newline]);
fprintf(fidout,['       EXTRAPOLATION_ORDER  ' num2str(EXOR) newline]);
fprintf(fidout,['       METHOD  GPW' newline]);
fprintf(fidout,['     &END QS' newline]);
fprintf(fidout,['     &MGRID' newline]);
fprintf(fidout,['       NGRIDS  4' newline]);
fprintf(fidout,['       CUTOFF     3.0000000000000000E+02' newline]);
fprintf(fidout,['       &RS_GRID' newline]);
fprintf(fidout,['         DISTRIBUTION_TYPE  REPLICATED' newline]);
fprintf(fidout,['       &END RS_GRID' newline]);
fprintf(fidout,['     &END MGRID' newline]);
fprintf(fidout,['     &XC' newline]);
fprintf(fidout,['       DENSITY_CUTOFF     1.0000000000000000E-10' newline]);
fprintf(fidout,['       GRADIENT_CUTOFF     1.0000000000000000E-10' newline]);
fprintf(fidout,['       TAU_CUTOFF     1.0000000000000000E-10' newline]);
fprintf(fidout,['       &XC_FUNCTIONAL  NO_SHORTCUT' newline]);
fprintf(fidout,['         &PBE  T' newline]);
fprintf(fidout,['         &END PBE' newline]);
fprintf(fidout,['       &END XC_FUNCTIONAL' newline]);
fprintf(fidout,['       &VDW_POTENTIAL' newline]);
fprintf(fidout,['         POTENTIAL_TYPE  PAIR_POTENTIAL' newline]);
fprintf(fidout,['         &PAIR_POTENTIAL' newline]);
fprintf(fidout,['           R_CUTOFF     8.0000000000000000E+00' newline]);
fprintf(fidout,['           TYPE  DFTD3' newline]);
fprintf(fidout,['           PARAMETER_FILE_NAME ./dftd3.dat' newline]);
fprintf(fidout,['           REFERENCE_FUNCTIONAL PBE' newline]);
fprintf(fidout,['           EPS_CN     1.0000000000000000E-02' newline]);
fprintf(fidout,['           CALCULATE_C9_TERM  T' newline]);
fprintf(fidout,['           REFERENCE_C9_TERM  T' newline]);
fprintf(fidout,['           LONG_RANGE_CORRECTION  T' newline]);
fprintf(fidout,['         &END PAIR_POTENTIAL' newline]);
fprintf(fidout,['       &END VDW_POTENTIAL' newline]);
fprintf(fidout,['     &END XC' newline]);
fprintf(fidout,['     &REAL_TIME_PROPAGATION' newline]);
fprintf(fidout,['       INITIAL_WFN  SCF_WFN' newline]);
fprintf(fidout,['     &END REAL_TIME_PROPAGATION' newline]);
fprintf(fidout,['   &END DFT' newline]);
fprintf(fidout,['   &SUBSYS' newline]);
fprintf(fidout,['     &CELL' newline]);
fprintf(fidout,[ABCOut newline]);
fprintf(fidout,['     &END CELL' newline]);
fprintf(fidout,['     &COORD' newline]);
for i = 1:length(CoordOut)
    fprintf(fidout,[CoordOut(i,:) newline]);
end
fprintf(fidout,['     &END COORD' newline]);
fprintf(fidout,['     &VELOCITY' newline]);
for i = 1:length(VelOut)
    fprintf(fidout,[VelOut(i,:) newline]);
end
fprintf(fidout,['     &END VELOCITY' newline]);
fprintf(fidout,['     &KIND Pt' newline]);
fprintf(fidout,['       BASIS_SET TZV-GTH-LDA-q18-very-confined' newline]);
fprintf(fidout,['       ELEMENT Pt' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q18' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND Ptb' newline]);
fprintf(fidout,['       BASIS_SET TZV-GTH-LDA-q18-very-confined' newline]);
fprintf(fidout,['       ELEMENT Pt' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q18' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND Pts' newline]);
fprintf(fidout,['       BASIS_SET TZV-GTH-LDA-q18-very-confined' newline]);
fprintf(fidout,['       ELEMENT Pt' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q18' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND Ptss' newline]);
fprintf(fidout,['       BASIS_SET TZV-GTH-LDA-q18-very-confined' newline]);
fprintf(fidout,['       ELEMENT Pt' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q18' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND H' newline]);
fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q1' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND O' newline]);
fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q6' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND Na' newline]);
fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q9' newline]);
fprintf(fidout,['     &END KIND' newline]);
fprintf(fidout,['     &KIND Cl' newline]);
fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
fprintf(fidout,['       POTENTIAL GTH-PBE-q7' newline]);
fprintf(fidout,['     &END KIND' newline]);
% fprintf(fidout,['     &TOPOLOGY' newline]);
% fprintf(fidout,['      COORD_FILE_NAME ' OutFldr '.xyz' newline]);
% fprintf(fidout,['      COORDINATE XYZ' newline]);
% fprintf(fidout,['     &END TOPOLOGY' newline]);
fprintf(fidout,['   &END SUBSYS' newline]);
fprintf(fidout,[' &END FORCE_EVAL' newline]);
fprintf(fidout,['' newline]);

fclose(fidout);