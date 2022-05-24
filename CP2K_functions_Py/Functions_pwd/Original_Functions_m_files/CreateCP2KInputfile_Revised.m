function CreateCP2KInputfile(Base, OutFldr, Vector, Wall, FixList, IonicSpecs, Calc, Velocities, Elems)

% write submission script

fidout = fopen([Base OutFldr '\' OutFldr '.inp'],'w');

fprintf(fidout,['# CP2K input file created by Matt Darby, Imperial College London for system ' OutFldr newline]);
fprintf(fidout,[' &GLOBAL' newline]);
fprintf(fidout,['   PRINT_LEVEL  LOW' newline]);
fprintf(fidout,['   PROJECT_NAME ' OutFldr newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'CP')
    fprintf(fidout,['   RUN_TYPE  MD' newline]);
elseif strcmp(Calc, 'GEO')
    fprintf(fidout,['   RUN_TYPE  GEO_OPT' newline]);
end
fprintf(fidout,['   WALLTIME ' num2str(Wall*60*60-600, '%i') newline]);
fprintf(fidout,['#   WALLTIME 1600' newline]);
fprintf(fidout,[' &END GLOBAL' newline]);
fprintf(fidout,[' &MOTION' newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'CP')
    fprintf(fidout,['   &MD' newline]);
    if strcmp(Calc, 'BOMD')
        fprintf(fidout,['     ENSEMBLE  NVT' newline]);
        fprintf(fidout,['	 &THERMOSTAT' newline]);
        fprintf(fidout,['       &CSVR' newline]);
        fprintf(fidout,['         TIMECON     1000' newline]);
        fprintf(fidout,['       &END CSVR' newline]);
        fprintf(fidout,['     &END THERMOSTAT' newline]);
    elseif strcmp(Calc, 'CP')
        fprintf(fidout,['     ENSEMBLE  LANGEVIN' newline]);
        fprintf(fidout,['     &LANGEVIN' newline]);
        x = [0.01; 0; 0; 3.75e-4];
        %         x = input('Enter array of STEPSIZE; EXOR; GAMMA; NOISY_GAMMA; ');
        fprintf(fidout,['       GAMMA     ' num2str(x(3)) newline]);
        fprintf(fidout,['       NOISY_GAMMA     ' num2str(x(4)) newline]);
        fprintf(fidout,['     &END LANGEVIN' newline]);
    end
    fprintf(fidout,['     STEPS  100000' newline]);
    fprintf(fidout,['     TIMESTEP     0.5' newline]);
    fprintf(fidout,['     TEMPERATURE     3.4000000000000000E+02' newline]);
    fprintf(fidout,['     TEMP_TOL     50.0000000000000000E+00' newline]);
    fprintf(fidout,['    &END MD' newline]);
elseif strcmp(Calc, 'GEO')
    fprintf(fidout,['   &GEO_OPT' newline]);
    fprintf(fidout,['     TYPE MINIMIZATION' newline]);
    fprintf(fidout,['     MAX_DR    1.0E-03' newline]);
    fprintf(fidout,['     MAX_FORCE 1.0E-03' newline]);
    fprintf(fidout,['     RMS_DR    1.0E-03' newline]);
    fprintf(fidout,['     RMS_FORCE 1.0E-03' newline]);
    fprintf(fidout,['     MAX_ITER 800' newline]);
    fprintf(fidout,['     OPTIMIZER BFGS' newline]);
    fprintf(fidout,['   &END GEO_OPT' newline]);
end
fprintf(fidout,['   &CONSTRAINT' newline]);
fprintf(fidout,['     &FIXED_ATOMS' newline]);
fprintf(fidout,['       LIST  ' num2str(FixList') newline]);
fprintf(fidout,['     &END FIXED_ATOMS' newline]);
fprintf(fidout,['   &END CONSTRAINT' newline]);
fprintf(fidout,['   &PRINT' newline]);
fprintf(fidout,['#     &VELOCITIES  ON' newline]);
fprintf(fidout,['#     &END VELOCITIES' newline]);
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
fprintf(fidout,['     &SCF' newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'GEO')
    fprintf(fidout,['       SCF_GUESS  RESTART' newline]);
    fprintf(fidout,['       MAX_SCF  800' newline]);
    fprintf(fidout,['       EPS_SCF     1.0E-06' newline]);
elseif strcmp(Calc, 'CP')
%     fprintf(fidout,['       SCF_GUESS  HISTORY_RESTART' newline]);
fprintf(fidout,['       SCF_GUESS  ATOMIC' newline]);
    fprintf(fidout,['       MAX_SCF_HISTORY  3' newline]);
    fprintf(fidout,['       EPS_SCF_HISTORY     1.0E-05' newline]);
    fprintf(fidout,['       MAX_SCF  800' newline]);
    fprintf(fidout,['       EPS_SCF     1.0E-06' newline]);
end
fprintf(fidout,['       EPS_DIIS     1.0E-01' newline]);
fprintf(fidout,['       &OT  T' newline]);
fprintf(fidout,['         MINIMIZER  DIIS' newline]);
fprintf(fidout,['         SAFE_DIIS  F' newline]);
fprintf(fidout,['         N_HISTORY_VEC  7' newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'GEO')
    fprintf(fidout,['         STEPSIZE     1.0000000000000001E-02' newline]);
elseif strcmp(Calc, 'CP')
    fprintf(fidout,['         STEPSIZE     ' num2str(x(1)) newline]);
end
fprintf(fidout,['         PRECONDITIONER  FULL_SINGLE_INVERSE' newline]);
fprintf(fidout,['         ENERGY_GAP     1.0000000000000000E-03' newline]);
fprintf(fidout,['       &END OT' newline]);
fprintf(fidout,['       &MIXING  T' newline]);
fprintf(fidout,['         METHOD  DIRECT_P_MIXING' newline]);
fprintf(fidout,['         ALPHA     3.000E-01' newline]);
fprintf(fidout,['       &END MIXING' newline]);
fprintf(fidout,['       &PRINT' newline]);
fprintf(fidout,['         &RESTART_HISTORY  SILENT' newline]);
if strcmp(Calc, 'BOMD')
    fprintf(fidout,['           BACKUP_COPIES  3' newline]);
elseif strcmp(Calc, 'CP')
    fprintf(fidout,['           BACKUP_COPIES  ' num2str(x(2)+2, '%i') newline]);
end
fprintf(fidout,['           &EACH' newline]);
fprintf(fidout,['           &END EACH' newline]);
fprintf(fidout,['         &END RESTART_HISTORY' newline]);
fprintf(fidout,['       &END PRINT' newline]);
fprintf(fidout,['     &END SCF' newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'CP')
    fprintf(fidout,['     &PRINT' newline]);
    fprintf(fidout,['         #&PDOS' newline]);
    fprintf(fidout,['         #  NLUMO 1000' newline]);
    fprintf(fidout,['         #  &EACH' newline]);
    fprintf(fidout,['         #    MD 500' newline]);
    fprintf(fidout,['         #  &END EACH' newline]);
    fprintf(fidout,['         #&END PDOS' newline]);
    fprintf(fidout,['         &V_HARTREE_CUBE' newline]);
    fprintf(fidout,['           STRIDE 1' newline]);
    fprintf(fidout,['           &EACH' newline]);
    fprintf(fidout,['             MD 500' newline]);
    fprintf(fidout,['           &END EACH' newline]);
    fprintf(fidout,['         &END V_HARTREE_CUBE' newline]);
    fprintf(fidout,['         &TOT_DENSITY_CUBE' newline]);
    fprintf(fidout,['           STRIDE 1' newline]);
    fprintf(fidout,['           &EACH' newline]);
    fprintf(fidout,['             MD 500' newline]);
    fprintf(fidout,['           &END EACH' newline]);
    fprintf(fidout,['         &END TOT_DENSITY_CUBE' newline]);
    fprintf(fidout,['         &E_DENSITY_CUBE' newline]);
    fprintf(fidout,['           STRIDE 1' newline]);
    fprintf(fidout,['           &EACH' newline]);
    fprintf(fidout,['             MD 500' newline]);
    fprintf(fidout,['           &END EACH' newline]);
    fprintf(fidout,['         &END E_DENSITY_CUBE' newline]);
    fprintf(fidout,['     &END PRINT' newline]);
end
fprintf(fidout,['     &QS' newline]);
if strcmp(Calc, 'BOMD') | strcmp(Calc, 'GEO')
    fprintf(fidout,['       EPS_DEFAULT     1.0E-11' newline]);
    fprintf(fidout,['       EXTRAPOLATION  ASPC' newline]);
    fprintf(fidout,['       EXTRAPOLATION_ORDER  0' newline]);
elseif strcmp(Calc, 'CP')
    fprintf(fidout,['       EPS_DEFAULT     1.0E-10' newline]);
    fprintf(fidout,['       EXTRAPOLATION  ASPC' newline]);
    fprintf(fidout,['       EXTRAPOLATION_ORDER  ' num2str(x(2), '%i') newline]);
end
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
fprintf(fidout,['       ABC     ' pad(num2str(Vector(1,1), '%.10f'),20) pad(num2str(Vector(2,2), '%.10f'),20) pad(num2str(Vector(3,3), '%.10f'),20) newline]);
fprintf(fidout,['     &END CELL' newline]);
if ~isempty(Velocities)
    fprintf(fidout,['     &VELOCITY' newline]);
    for i = 1:size(Velocities,1)
        fprintf(fidout,['       ' pad(num2str(Velocities(i,1), '%.10f'),20) pad(num2str(Velocities(i,2), '%.10f'),20) pad(num2str(Velocities(i,3), '%.10f'),20)  newline]);
    end
    fprintf(fidout,['     &END VELOCITY' newline]);
end
for bp = 1:size(Elems{1},1)
    fprintf(fidout,[strjoin(['     &KIND ' Elems{2}(bp,:)]) newline]);
    if strcmp(Elems{1}(bp,:), 'Pt')
        fprintf(fidout,['       ELEMENT Pt' newline]);
        fprintf(fidout,['       BASIS_SET TZV-GTH-LDA-q18-very-confined' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q18' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Al')
        fprintf(fidout,['       ELEMENT Al' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-MOLOPT-SR-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q3' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Ag')
        fprintf(fidout,['       ELEMENT Ag' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-MOLOPT-SR-GTH-q11' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q11' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Au')
        fprintf(fidout,['       ELEMENT Au' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-MOLOPT-SR-GTH-q11' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q11' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Ir')
        fprintf(fidout,['       ELEMENT Ir' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-MOLOPT-SR-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q17' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Ru')
        fprintf(fidout,['       ELEMENT Ru' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-MOLOPT-SR-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q16' newline]);
    elseif strcmp(Elems{1}(bp,:), 'F')
        fprintf(fidout,['       ELEMENT F' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q7' newline]);
    elseif strcmp(Elems{1}(bp,:), 'O')
        fprintf(fidout,['       ELEMENT O' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q6' newline]);
    elseif strcmp(Elems{1}(bp,:), 'H')
        fprintf(fidout,['       ELEMENT H' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q1' newline]);
    elseif strcmp(Elems{1}(bp,:), 'Cl')
        fprintf(fidout,['       ELEMENT Cl' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q7' newline]);
        elseif strcmp(Elems{1}(bp,:), 'Mo')
        fprintf(fidout,['       ELEMENT Mo' newline]);
        fprintf(fidout,['       BASIS_SET TZVP-GTH' newline]);
        fprintf(fidout,['       POTENTIAL GTH-PBE-q7' newline]);
    end
    
    fprintf(fidout,['     &END KIND' newline]);
end
fprintf(fidout,['     &TOPOLOGY' newline]);
fprintf(fidout,['      COORD_FILE_NAME ' OutFldr '.xyz' newline]);
fprintf(fidout,['      COORDINATE XYZ' newline]);
fprintf(fidout,['     &END TOPOLOGY' newline]);
fprintf(fidout,['   &END SUBSYS' newline]);
fprintf(fidout,[' &END FORCE_EVAL' newline]);
fprintf(fidout,['' newline]);

fclose(fidout);

copyfile([Base OutFldr '\' OutFldr '.inp'], [Base OutFldr '\' OutFldr '-1.restart'])