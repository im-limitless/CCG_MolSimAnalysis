import numpy as np
    
def CreateCP2KInputfile(Base = None,OutFldr = None,Vector = None,Wall = None,FixList = None,IonicSpecs = None,Calc = None,Velocities = None,Elems = None): 
    # write submission script
    
    fidout = open(np.array([Base,OutFldr,'\',OutFldr,'.inp']),'w')
    fidout.write(np.array(['# CP2K input file created by Matt Darby, Imperial College London for system ',OutFldr,newline]) % ())
    fidout.write(np.array([' &GLOBAL',newline]) % ())
    fidout.write(np.array(['   PRINT_LEVEL  LOW',newline]) % ())
    fidout.write(np.array(['   PROJECT_NAME ',OutFldr,newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('CP')):
        fidout.write(np.array(['   RUN_TYPE  MD',newline]) % ())
    else:
        if str(Calc) == str('GEO'):
            fidout.write(np.array(['   RUN_TYPE  GEO_OPT',newline]) % ())
    
    fidout.write(np.array(['   WALLTIME ',num2str(Wall * 60 * 60 - 600,'%i'),newline]) % ())
    fidout.write(np.array(['#   WALLTIME 1600',newline]) % ())
    fidout.write(np.array([' &END GLOBAL',newline]) % ())
    fidout.write(np.array([' &MOTION',newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('CP')):
        fidout.write(np.array(['   &MD',newline]) % ())
        if str(Calc) == str('BOMD'):
            fidout.write(np.array(['     ENSEMBLE  NVT',newline]) % ())
            fidout.write(np.array(['	 &THERMOSTAT',newline]) % ())
            fidout.write(np.array(['       &CSVR',newline]) % ())
            fidout.write(np.array(['         TIMECON     1000',newline]) % ())
            fidout.write(np.array(['       &END CSVR',newline]) % ())
            fidout.write(np.array(['     &END THERMOSTAT',newline]) % ())
        else:
            if str(Calc) == str('CP'):
                fidout.write(np.array(['     ENSEMBLE  LANGEVIN',newline]) % ())
                fidout.write(np.array(['     &LANGEVIN',newline]) % ())
                x = np.array([[0.01],[0],[0],[0.000375]])
                #         x = input('Enter array of STEPSIZE; EXOR; GAMMA; NOISY_GAMMA; ');
                fidout.write(np.array(['       GAMMA     ',num2str(x(3)),newline]) % ())
                fidout.write(np.array(['       NOISY_GAMMA     ',num2str(x(4)),newline]) % ())
                fidout.write(np.array(['     &END LANGEVIN',newline]) % ())
        fidout.write(np.array(['     STEPS  100000',newline]) % ())
        fidout.write(np.array(['     TIMESTEP     0.5',newline]) % ())
        fidout.write(np.array(['     TEMPERATURE     3.4000000000000000E+02',newline]) % ())
        fidout.write(np.array(['     TEMP_TOL     50.0000000000000000E+00',newline]) % ())
        fidout.write(np.array(['    &END MD',newline]) % ())
    else:
        if str(Calc) == str('GEO'):
            fidout.write(np.array(['   &GEO_OPT',newline]) % ())
            fidout.write(np.array(['     TYPE MINIMIZATION',newline]) % ())
            fidout.write(np.array(['     MAX_DR    1.0E-03',newline]) % ())
            fidout.write(np.array(['     MAX_FORCE 1.0E-03',newline]) % ())
            fidout.write(np.array(['     RMS_DR    1.0E-03',newline]) % ())
            fidout.write(np.array(['     RMS_FORCE 1.0E-03',newline]) % ())
            fidout.write(np.array(['     MAX_ITER 800',newline]) % ())
            fidout.write(np.array(['     OPTIMIZER BFGS',newline]) % ())
            fidout.write(np.array(['   &END GEO_OPT',newline]) % ())
    
    fidout.write(np.array(['   &CONSTRAINT',newline]) % ())
    fidout.write(np.array(['     &FIXED_ATOMS',newline]) % ())
    fidout.write(np.array(['       LIST  ',num2str(np.transpose(FixList)),newline]) % ())
    fidout.write(np.array(['     &END FIXED_ATOMS',newline]) % ())
    fidout.write(np.array(['   &END CONSTRAINT',newline]) % ())
    fidout.write(np.array(['   &PRINT',newline]) % ())
    fidout.write(np.array(['#     &VELOCITIES  ON',newline]) % ())
    fidout.write(np.array(['#     &END VELOCITIES',newline]) % ())
    fidout.write(np.array(['     &RESTART  SILENT',newline]) % ())
    fidout.write(np.array(['       ADD_LAST  NUMERIC',newline]) % ())
    fidout.write(np.array(['       &EACH',newline]) % ())
    fidout.write(np.array(['         MD  1',newline]) % ())
    fidout.write(np.array(['       &END EACH',newline]) % ())
    fidout.write(np.array(['     &END RESTART',newline]) % ())
    fidout.write(np.array(['   &END PRINT',newline]) % ())
    fidout.write(np.array([' &END MOTION',newline]) % ())
    fidout.write(np.array([' &FORCE_EVAL',newline]) % ())
    fidout.write(np.array(['   METHOD  QS',newline]) % ())
    #     fprintf(fidout,['   &EXTERNAL_POTENTIAL' newline]);
#     fprintf(fidout,['    ATOMS_LIST ' num2str(IonicSpecs') newline]);
#     fprintf(fidout,['    FUNCTION (1.0E-6)*((Z-' num2str(Vector(3,3)*1.88973/2) ')^4)' newline]);
#     fprintf(fidout,['   &END EXTERNAL_POTENTIAL' newline]);
    fidout.write(np.array(['   &DFT',newline]) % ())
    fidout.write(np.array(['     BASIS_SET_FILE_NAME ./GTH_BASIS_SETS',newline]) % ())
    fidout.write(np.array(['     POTENTIAL_FILE_NAME ./GTH_POTENTIALS',newline]) % ())
    fidout.write(np.array(['     &SCF',newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('GEO')):
        fidout.write(np.array(['       SCF_GUESS  RESTART',newline]) % ())
        fidout.write(np.array(['       MAX_SCF  800',newline]) % ())
        fidout.write(np.array(['       EPS_SCF     1.0E-06',newline]) % ())
    else:
        if str(Calc) == str('CP'):
            #     fprintf(fidout,['       SCF_GUESS  HISTORY_RESTART' newline]);
            fidout.write(np.array(['       SCF_GUESS  ATOMIC',newline]) % ())
            fidout.write(np.array(['       MAX_SCF_HISTORY  3',newline]) % ())
            fidout.write(np.array(['       EPS_SCF_HISTORY     1.0E-05',newline]) % ())
            fidout.write(np.array(['       MAX_SCF  800',newline]) % ())
            fidout.write(np.array(['       EPS_SCF     1.0E-06',newline]) % ())
    
    fidout.write(np.array(['       EPS_DIIS     1.0E-01',newline]) % ())
    fidout.write(np.array(['       &OT  T',newline]) % ())
    fidout.write(np.array(['         MINIMIZER  DIIS',newline]) % ())
    fidout.write(np.array(['         SAFE_DIIS  F',newline]) % ())
    fidout.write(np.array(['         N_HISTORY_VEC  7',newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('GEO')):
        fidout.write(np.array(['         STEPSIZE     1.0000000000000001E-02',newline]) % ())
    else:
        if str(Calc) == str('CP'):
            fidout.write(np.array(['         STEPSIZE     ',num2str(x(1)),newline]) % ())
    
    fidout.write(np.array(['         PRECONDITIONER  FULL_SINGLE_INVERSE',newline]) % ())
    fidout.write(np.array(['         ENERGY_GAP     1.0000000000000000E-03',newline]) % ())
    fidout.write(np.array(['       &END OT',newline]) % ())
    fidout.write(np.array(['       &MIXING  T',newline]) % ())
    fidout.write(np.array(['         METHOD  DIRECT_P_MIXING',newline]) % ())
    fidout.write(np.array(['         ALPHA     3.000E-01',newline]) % ())
    fidout.write(np.array(['       &END MIXING',newline]) % ())
    fidout.write(np.array(['       &PRINT',newline]) % ())
    fidout.write(np.array(['         &RESTART_HISTORY  SILENT',newline]) % ())
    if str(Calc) == str('BOMD'):
        fidout.write(np.array(['           BACKUP_COPIES  3',newline]) % ())
    else:
        if str(Calc) == str('CP'):
            fidout.write(np.array(['           BACKUP_COPIES  ',num2str(x(2) + 2,'%i'),newline]) % ())
    
    fidout.write(np.array(['           &EACH',newline]) % ())
    fidout.write(np.array(['           &END EACH',newline]) % ())
    fidout.write(np.array(['         &END RESTART_HISTORY',newline]) % ())
    fidout.write(np.array(['       &END PRINT',newline]) % ())
    fidout.write(np.array(['     &END SCF',newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('CP')):
        fidout.write(np.array(['     &PRINT',newline]) % ())
        fidout.write(np.array(['         #&PDOS',newline]) % ())
        fidout.write(np.array(['         #  NLUMO 1000',newline]) % ())
        fidout.write(np.array(['         #  &EACH',newline]) % ())
        fidout.write(np.array(['         #    MD 500',newline]) % ())
        fidout.write(np.array(['         #  &END EACH',newline]) % ())
        fidout.write(np.array(['         #&END PDOS',newline]) % ())
        fidout.write(np.array(['         &V_HARTREE_CUBE',newline]) % ())
        fidout.write(np.array(['           STRIDE 1',newline]) % ())
        fidout.write(np.array(['           &EACH',newline]) % ())
        fidout.write(np.array(['             MD 500',newline]) % ())
        fidout.write(np.array(['           &END EACH',newline]) % ())
        fidout.write(np.array(['         &END V_HARTREE_CUBE',newline]) % ())
        fidout.write(np.array(['         &TOT_DENSITY_CUBE',newline]) % ())
        fidout.write(np.array(['           STRIDE 1',newline]) % ())
        fidout.write(np.array(['           &EACH',newline]) % ())
        fidout.write(np.array(['             MD 500',newline]) % ())
        fidout.write(np.array(['           &END EACH',newline]) % ())
        fidout.write(np.array(['         &END TOT_DENSITY_CUBE',newline]) % ())
        fidout.write(np.array(['         &E_DENSITY_CUBE',newline]) % ())
        fidout.write(np.array(['           STRIDE 1',newline]) % ())
        fidout.write(np.array(['           &EACH',newline]) % ())
        fidout.write(np.array(['             MD 500',newline]) % ())
        fidout.write(np.array(['           &END EACH',newline]) % ())
        fidout.write(np.array(['         &END E_DENSITY_CUBE',newline]) % ())
        fidout.write(np.array(['     &END PRINT',newline]) % ())
    
    fidout.write(np.array(['     &QS',newline]) % ())
    if np.logical_or(str(Calc) == str('BOMD'),str(Calc) == str('GEO')):
        fidout.write(np.array(['       EPS_DEFAULT     1.0E-11',newline]) % ())
        fidout.write(np.array(['       EXTRAPOLATION  ASPC',newline]) % ())
        fidout.write(np.array(['       EXTRAPOLATION_ORDER  0',newline]) % ())
    else:
        if str(Calc) == str('CP'):
            fidout.write(np.array(['       EPS_DEFAULT     1.0E-10',newline]) % ())
            fidout.write(np.array(['       EXTRAPOLATION  ASPC',newline]) % ())
            fidout.write(np.array(['       EXTRAPOLATION_ORDER  ',num2str(x(2),'%i'),newline]) % ())
    
    fidout.write(np.array(['       METHOD  GPW',newline]) % ())
    fidout.write(np.array(['     &END QS',newline]) % ())
    fidout.write(np.array(['     &MGRID',newline]) % ())
    fidout.write(np.array(['       NGRIDS  4',newline]) % ())
    fidout.write(np.array(['       CUTOFF     3.0000000000000000E+02',newline]) % ())
    fidout.write(np.array(['       &RS_GRID',newline]) % ())
    fidout.write(np.array(['         DISTRIBUTION_TYPE  REPLICATED',newline]) % ())
    fidout.write(np.array(['       &END RS_GRID',newline]) % ())
    fidout.write(np.array(['     &END MGRID',newline]) % ())
    fidout.write(np.array(['     &XC',newline]) % ())
    fidout.write(np.array(['       DENSITY_CUTOFF     1.0000000000000000E-10',newline]) % ())
    fidout.write(np.array(['       GRADIENT_CUTOFF     1.0000000000000000E-10',newline]) % ())
    fidout.write(np.array(['       TAU_CUTOFF     1.0000000000000000E-10',newline]) % ())
    fidout.write(np.array(['       &XC_FUNCTIONAL  NO_SHORTCUT',newline]) % ())
    fidout.write(np.array(['         &PBE  T',newline]) % ())
    fidout.write(np.array(['         &END PBE',newline]) % ())
    fidout.write(np.array(['       &END XC_FUNCTIONAL',newline]) % ())
    fidout.write(np.array(['       &VDW_POTENTIAL',newline]) % ())
    fidout.write(np.array(['         POTENTIAL_TYPE  PAIR_POTENTIAL',newline]) % ())
    fidout.write(np.array(['         &PAIR_POTENTIAL',newline]) % ())
    fidout.write(np.array(['           R_CUTOFF     8.0000000000000000E+00',newline]) % ())
    fidout.write(np.array(['           TYPE  DFTD3',newline]) % ())
    fidout.write(np.array(['           PARAMETER_FILE_NAME ./dftd3.dat',newline]) % ())
    fidout.write(np.array(['           REFERENCE_FUNCTIONAL PBE',newline]) % ())
    fidout.write(np.array(['           EPS_CN     1.0000000000000000E-02',newline]) % ())
    fidout.write(np.array(['           CALCULATE_C9_TERM  T',newline]) % ())
    fidout.write(np.array(['           REFERENCE_C9_TERM  T',newline]) % ())
    fidout.write(np.array(['           LONG_RANGE_CORRECTION  T',newline]) % ())
    fidout.write(np.array(['         &END PAIR_POTENTIAL',newline]) % ())
    fidout.write(np.array(['       &END VDW_POTENTIAL',newline]) % ())
    fidout.write(np.array(['     &END XC',newline]) % ())
    fidout.write(np.array(['     &REAL_TIME_PROPAGATION',newline]) % ())
    fidout.write(np.array(['       INITIAL_WFN  SCF_WFN',newline]) % ())
    fidout.write(np.array(['     &END REAL_TIME_PROPAGATION',newline]) % ())
    fidout.write(np.array(['   &END DFT',newline]) % ())
    fidout.write(np.array(['   &SUBSYS',newline]) % ())
    fidout.write(np.array(['     &CELL',newline]) % ())
    fidout.write(np.array(['       ABC     ',pad(num2str(Vector(1,1),'%.10f'),20),pad(num2str(Vector(2,2),'%.10f'),20),pad(num2str(Vector(3,3),'%.10f'),20),newline]) % ())
    fidout.write(np.array(['     &END CELL',newline]) % ())
    if not len(Velocities)==0 :
        fidout.write(np.array(['     &VELOCITY',newline]) % ())
        for i in np.arange(1,Velocities.shape[1-1]+1).reshape(-1):
            fidout.write(np.array(['       ',pad(num2str(Velocities(i,1),'%.10f'),20),pad(num2str(Velocities(i,2),'%.10f'),20),pad(num2str(Velocities(i,3),'%.10f'),20),newline]) % ())
        fidout.write(np.array(['     &END VELOCITY',newline]) % ())
    
    for bp in np.arange(1,Elems[0].shape[1-1]+1).reshape(-1):
        fidout.write(np.array([strjoin(np.array(['     &KIND ',Elems[2](bp,:)])),newline]) % ())
        if str(Elems[0](bp,:)) == str('Pt'):
            fidout.write(np.array(['       ELEMENT Pt',newline]) % ())
            fidout.write(np.array(['       BASIS_SET TZV-GTH-LDA-q18-very-confined',newline]) % ())
            fidout.write(np.array(['       POTENTIAL GTH-PBE-q18',newline]) % ())
        else:
            if str(Elems[0](bp,:)) == str('Al'):
                fidout.write(np.array(['       ELEMENT Al',newline]) % ())
                fidout.write(np.array(['       BASIS_SET TZVP-MOLOPT-SR-GTH',newline]) % ())
                fidout.write(np.array(['       POTENTIAL GTH-PBE-q3',newline]) % ())
            else:
                if str(Elems[0](bp,:)) == str('Ag'):
                    fidout.write(np.array(['       ELEMENT Ag',newline]) % ())
                    fidout.write(np.array(['       BASIS_SET TZVP-MOLOPT-SR-GTH-q11',newline]) % ())
                    fidout.write(np.array(['       POTENTIAL GTH-PBE-q11',newline]) % ())
                else:
                    if str(Elems[0](bp,:)) == str('Au'):
                        fidout.write(np.array(['       ELEMENT Au',newline]) % ())
                        fidout.write(np.array(['       BASIS_SET TZVP-MOLOPT-SR-GTH-q11',newline]) % ())
                        fidout.write(np.array(['       POTENTIAL GTH-PBE-q11',newline]) % ())
                    else:
                        if str(Elems[0](bp,:)) == str('Ir'):
                            fidout.write(np.array(['       ELEMENT Ir',newline]) % ())
                            fidout.write(np.array(['       BASIS_SET TZVP-MOLOPT-SR-GTH',newline]) % ())
                            fidout.write(np.array(['       POTENTIAL GTH-PBE-q17',newline]) % ())
                        else:
                            if str(Elems[0](bp,:)) == str('Ru'):
                                fidout.write(np.array(['       ELEMENT Ru',newline]) % ())
                                fidout.write(np.array(['       BASIS_SET TZVP-MOLOPT-SR-GTH',newline]) % ())
                                fidout.write(np.array(['       POTENTIAL GTH-PBE-q16',newline]) % ())
                            else:
                                if str(Elems[0](bp,:)) == str('F'):
                                    fidout.write(np.array(['       ELEMENT F',newline]) % ())
                                    fidout.write(np.array(['       BASIS_SET TZVP-GTH',newline]) % ())
                                    fidout.write(np.array(['       POTENTIAL GTH-PBE-q7',newline]) % ())
                                else:
                                    if str(Elems[0](bp,:)) == str('O'):
                                        fidout.write(np.array(['       ELEMENT O',newline]) % ())
                                        fidout.write(np.array(['       BASIS_SET TZVP-GTH',newline]) % ())
                                        fidout.write(np.array(['       POTENTIAL GTH-PBE-q6',newline]) % ())
                                    else:
                                        if str(Elems[0](bp,:)) == str('H'):
                                            fidout.write(np.array(['       ELEMENT H',newline]) % ())
                                            fidout.write(np.array(['       BASIS_SET TZVP-GTH',newline]) % ())
                                            fidout.write(np.array(['       POTENTIAL GTH-PBE-q1',newline]) % ())
                                        else:
                                            if str(Elems[0](bp,:)) == str('Cl'):
                                                fidout.write(np.array(['       ELEMENT Cl',newline]) % ())
                                                fidout.write(np.array(['       BASIS_SET TZVP-GTH',newline]) % ())
                                                fidout.write(np.array(['       POTENTIAL GTH-PBE-q7',newline]) % ())
                                            else:
                                                if str(Elems[0](bp,:)) == str('Mo'):
                                                    fidout.write(np.array(['       ELEMENT Mo',newline]) % ())
                                                    fidout.write(np.array(['       BASIS_SET TZVP-GTH',newline]) % ())
                                                    fidout.write(np.array(['       POTENTIAL GTH-PBE-q7',newline]) % ())
        fidout.write(np.array(['     &END KIND',newline]) % ())
    
    fidout.write(np.array(['     &TOPOLOGY',newline]) % ())
    fidout.write(np.array(['      COORD_FILE_NAME ',OutFldr,'.xyz',newline]) % ())
    fidout.write(np.array(['      COORDINATE XYZ',newline]) % ())
    fidout.write(np.array(['     &END TOPOLOGY',newline]) % ())
    fidout.write(np.array(['   &END SUBSYS',newline]) % ())
    fidout.write(np.array([' &END FORCE_EVAL',newline]) % ())
    fidout.write(np.array(['',newline]) % ())
    fidout.close()
    copyfile(np.array([Base,OutFldr,'\',OutFldr,'.inp']),np.array([Base,OutFldr,'\',OutFldr,'-1.restart']))