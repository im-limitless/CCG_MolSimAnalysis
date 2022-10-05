function CreateSubmissionFile(Cluster, Base, OutFldr, Time, CPUs, CPUpn, VASPsol)

% write submission script

    if strcmp('Legion', Cluster) | strcmp('Grace', Cluster) | strcmp('Thomas', Cluster)
        fidout = fopen([Base OutFldr '\vasp_parallel.qs'],'w');
        
        fprintf(fidout,['#!/bin/bash -l' newline]);
        fprintf(fidout,['#$ -S /bin/bash' newline]);
        fprintf(fidout,['# Set wallclock time (format hours:minutes:seconds).' newline]);
        %     fprintf(fidout,['#$ -l paid=1' newline]);
        if strcmp('Legion', Cluster)
            fprintf(fidout,['#$ -l h_rt=' num2str(Time) ':00:00' newline]);
        else if strcmp('Grace', Cluster) | strcmp('Thomas', Cluster)
                fprintf(fidout,['#$ -l h_rt=' num2str(Time) ':00:00' newline]);
            end
        end
        fprintf(fidout,['# Set the name of the job.' newline]);
        fprintf(fidout,['#$ -N ' OutFldr newline]);
        fprintf(fidout,['#$ -M matthew.darby.13@ucl.ac.uk' newline]);
        fprintf(fidout,['#$ -m abe' newline]);
        fprintf(fidout,['# Select parallel environment.' newline]);
        fprintf(fidout,['#$ -pe qlc ' num2str(CPUs) newline]);
        fprintf(fidout,['# Set the working directory to the current one.' newline]);
        fprintf(fidout,['#$ -cwd ' newline]);
        fprintf(fidout,['' newline]);
        if strcmp('Legion', Cluster)
            if strcmp('No', VASPsol)
                fprintf(fidout,['APP=/shared/ucl/apps/vasp/5.4.1-p2-p2-vtst3-r160/intel-2015-update2/bin/vasp_std' newline]);
            elseif strcmp('Yes', VASPsol)
                fprintf(fidout,['module unload compilers mpi' newline]);
                fprintf(fidout,['module load compilers/intel/2017/update3' newline]);
                fprintf(fidout,['module load mpi/intel/2017/update3/intel' newline]);
                fprintf(fidout,['APP=~/VASP/vasp.5.4.1_VTSTr178_sol/bin/vasp_std' newline]);
            end
        elseif strcmp('Grace', Cluster)
            if strcmp('No', VASPsol)
                fprintf(fidout,['APP=~/VASP/vasp.5.4.1_VTST/bin/vasp_std' newline]);
            elseif strcmp('Yes', VASPsol)
                fprintf(fidout,['module unload compilers mpi' newline]);
                fprintf(fidout,['module load compilers/intel/2017/update3' newline]);
                fprintf(fidout,['module load mpi/intel/2017/update3/intel' newline]);
                fprintf(fidout,['APP=~/VASP/vasp.5.4.1_VTSTr178_sol/bin/vasp_std' newline]);
            end
        elseif strcmp('Thomas', Cluster)
                    fprintf(fidout,['#$ -P Gold' newline]);
                    fprintf(fidout,['#$ -A UCL_chemE' newline]);
                    fprintf(fidout,['APP=/home/UCL/apps/vasp/5.4.1-p2-p2-vtst3-r160/intel-2015-update2/bin/vasp_std' newline]);
        end
        fprintf(fidout,['WORKDIR=$SGE_O_WORKDIR' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# Count number of nodes from SGE_NODEFILE' newline]);
        fprintf(fidout,['NODECOUNT=`sort -u $SGE_NODEFILE | wc -l`' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# Count number of processors' newline]);
        fprintf(fidout,['PROCCOUNT=`cat $SGE_NODEFILE | wc -l`' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# execute Vasp' newline]);
        fprintf(fidout,['cd $WORKDIR' newline]);
        fprintf(fidout,['cat $PBS_NODEFILE > nodefile.out' newline]);
        fprintf(fidout,['gerun $APP > vasp.out' newline]);
        fprintf(fidout,['' newline]);
        
        fclose(fidout);
    end
    
    if strcmp('IBChem', Cluster)
        fidout = fopen([Base OutFldr '\vasp_parallel.qs'],'w');
        fprintf(fidout,['#!/bin/bash -l' newline]);
        fprintf(fidout,['#$ -S /bin/bash' newline]);
        fprintf(fidout,['# Set wallclock time (format hours:minutes:seconds).' newline]);
        fprintf(fidout,['#$ -l h_rt=' num2str(Time) ':00:00' newline]);
        fprintf(fidout,['# Set the name of the job.' newline]);
        fprintf(fidout,['#$ -N ' OutFldr newline]);
        fprintf(fidout,['# Select parallel environment.' newline]);
        fprintf(fidout,['#$ -pe openmpi ' num2str(CPUs) newline]);
        fprintf(fidout,['# Set the working directory to the current one.' newline]);
        fprintf(fidout,['#$ -cwd ' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['APP=/usr/local/vasp/vasp-5.4.1-p3-tst' newline]);
        fprintf(fidout,['WORKDIR=$SGE_O_WORKDIR' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# Count number of nodes from SGE_NODEFILE' newline]);
        fprintf(fidout,['NODECOUNT=`sort -u $SGE_NODEFILE | wc -l`' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# Count number of processors' newline]);
        fprintf(fidout,['PROCCOUNT=`cat $SGE_NODEFILE | wc -l`' newline]);
        fprintf(fidout,['' newline]);
        fprintf(fidout,['# execute Vasp' newline]);
        fprintf(fidout,['cd $WORKDIR' newline]);
        fprintf(fidout,['cat $PBS_NODEFILE > nodefile.out' newline]);
        fprintf(fidout,['mpirun $APP > vasp.out' newline]);
        fprintf(fidout,['' newline]);
        
        fclose(fidout);
    end
    if strcmp('LRZ', Cluster)
        fidout = fopen([Base OutFldr '\cp2k_parallel.qs'],'w');
        fprintf(fidout, '%s',['#!/bin/bash -l' newline]);
        fprintf(fidout, '%s',['# Job Name and Files (also --job-name) ' newline]);
        fprintf(fidout, '%s',['#SBATCH -J  ' OutFldr newline]);
        fprintf(fidout, '%s',['#Output and error (also --output, --error): ' newline]);
        fprintf(fidout, '%s',['#SBATCH -o ./%x.%j.out' newline]);
        fprintf(fidout, '%s',['#SBATCH -e ./%x.%j.err' newline]);
%         fprintf(fidout, '%s',['#SBATCH --ntasks=' num2str(CPUs) newline]);
        fprintf(fidout, '%s',['#SBATCH --nodes=' num2str(CPUs/48) newline]);
        fprintf(fidout, '%s',['#SBATCH --ntasks-per-node=' num2str(CPUpn) newline]);
        fprintf(fidout, '%s',['#Initial working directory (also --chdir):' newline]);
        fprintf(fidout, '%s',['#SBATCH -D ./' newline]);
        fprintf(fidout, '%s',['#Notification and type' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-type=BEGIN' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-type=END' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-user=m.darby@imperial.ac.uk' newline]);
        fprintf(fidout, '%s',['# Wall clock limit:' newline]);
        fprintf(fidout, '%s',['#SBATCH --time=' num2str(Time) ':00:00' newline]);
        fprintf(fidout, '%s',['#SBATCH --no-requeue' newline]);
        fprintf(fidout, '%s',['#Setup of execution environment' newline]);
        fprintf(fidout, '%s',['#SBATCH --export=NONE' newline]);
        fprintf(fidout, '%s',['#SBATCH --get-user-env' newline]);
        fprintf(fidout, '%s',['#SBATCH --account=pn29la' newline]);
        if CPUs < 816
            fprintf(fidout, '%s',['#SBATCH --partition=micro' newline]);
        elseif CPUs < 36912 &  CPUs >= 816
            fprintf(fidout, '%s',['#SBATCH --partition=general' newline]);
        elseif CPUs >= 36912
            fprintf(fidout, '%s',['#SBATCH --partition=large' newline]);
        end
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['module load slurm_setup' newline]);
        fprintf(fidout, '%s',['module load list' newline]);
        fprintf(fidout, '%s',['module load cp2k/7.1-gcc8-impi' newline]);
        fprintf(fidout, '%s',['cd ./' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['mpiexec -n ' num2str(CPUpn*CPUs/48) ' -ppn ' num2str(CPUpn) '  /dss/dsshome1/lrz/sys/spack/.tmp.test.molcas/master_010221/opt/skylake/cp2k/7.1-gcc-mmoxx2f/bin/cp2k.psmp  ' OutFldr '-1.restart 2>&1 > out.log' newline]);
%         fprintf(fidout, '%s',['mpiexec -n ' num2str(CPUs) ' cp2k.popt '  OutFldr '-1.restart >& out.log' newline]);
%         fprintf(fidout, '%s',['mpiexec -n ' num2str(CPUs) ' cp2k.popt '  OutFldr '.inp >& out.log' newline]);
        fprintf(fidout, '%s',['' newline]);
        
        fclose(fidout);
%     end
%     if strcmp('LRZ-test', Cluster)
        
        
        fidout = fopen([Base OutFldr '\test.qs'],'w');
        fprintf(fidout, '%s',['#!/bin/bash -l' newline]);
        fprintf(fidout, '%s',['# Job Name and Files (also --job-name) ' newline]);
        fprintf(fidout, '%s',['#SBATCH -J  ' OutFldr newline]);
        fprintf(fidout, '%s',['#Output and error (also --output, --error): ' newline]);
        fprintf(fidout, '%s', ['#SBATCH -o ./%x.%j.out' newline]);
        fprintf(fidout, '%s',['#SBATCH -e ./%x.%j.err' newline]);
%         fprintf(fidout, '%s',['#SBATCH --ntasks=768' newline]);
        fprintf(fidout, '%s',['#SBATCH --nodes=16' newline]);
        fprintf(fidout, '%s',['#SBATCH --ntasks-per-node=' num2str(CPUpn) newline]);
        fprintf(fidout, '%s',['#Initial working directory (also --chdir):' newline]);
        fprintf(fidout, '%s',['#SBATCH -D ./' newline]);
        fprintf(fidout, '%s',['#Notification and type' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-type=BEGIN' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-type=END' newline]);
%         fprintf(fidout, '%s',['#SBATCH --mail-user=m.darby@imperial.ac.uk' newline]);
        fprintf(fidout, '%s',['# Wall clock limit:' newline]);
        fprintf(fidout, '%s',['#SBATCH --time=00:30:00' newline]);
        fprintf(fidout, '%s',['#SBATCH --no-requeue' newline]);
        fprintf(fidout, '%s',['#Setup of execution environment' newline]);
        fprintf(fidout, '%s',['#SBATCH --export=NONE' newline]);
        fprintf(fidout, '%s',['#SBATCH --get-user-env' newline]);
        fprintf(fidout, '%s',['#SBATCH --account=pn29la' newline]);
        fprintf(fidout, '%s',['#SBATCH --partition=test' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['module load slurm_setup' newline]);
        fprintf(fidout, '%s',['module load list' newline]);
        fprintf(fidout, '%s',['module load cp2k/7.1-gcc8-impi' newline]);
%         fprintf(fidout, '%s',['module load cp2k' newline]);
        
        fprintf(fidout, '%s',['cd ./' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['cd ./' newline]);
%         fprintf(fidout, '%s',['mpiexec -n 768 cp2k.popt ' OutFldr '.inp >& out.log' newline]);
        fprintf(fidout, '%s',['mpiexec -n ' num2str(CPUpn*16) ' -ppn ' num2str(CPUpn) '  /dss/dsshome1/lrz/sys/spack/.tmp.test.molcas/master_010221/opt/skylake/cp2k/7.1-gcc-mmoxx2f/bin/cp2k.psmp  ' OutFldr '.inp 2>&1 > out.log' newline]);
        fprintf(fidout, '%s',['' newline]);
        
        fclose(fidout);
        
    end
    
    if strcmp('ARCHER', Cluster)
        fidout = fopen([Base OutFldr '\cp2k_parallel.qs'],'w');
        fprintf(fidout, '%s',['#!/bin/bash -l' newline]);
        fprintf(fidout, '%s',['# Job Name and Files (also --job-name) ' newline]);
        fprintf(fidout, '%s',['#SBATCH -J  ' OutFldr newline]);
        fprintf(fidout, '%s',['#Output and error (also --output, --error): ' newline]);
        fprintf(fidout, '%s',['#SBATCH -o ./%x.%j.out' newline]);
        fprintf(fidout, '%s',['#SBATCH -e ./%x.%j.err' newline]);
        fprintf(fidout, '%s',['#SBATCH --nodes=' num2str(CPUs/128) newline]);
        fprintf(fidout, '%s',['#SBATCH --ntasks-per-node=128' newline]);
        fprintf(fidout, '%s',['#SBATCH --cpus-per-task=1' newline]);
        fprintf(fidout, '%s',['#Initial working directory (also --chdir):' newline]);
        fprintf(fidout, '%s',['#SBATCH -D ./' newline]);
        fprintf(fidout, '%s',['# Wall clock limit:' newline]);
        fprintf(fidout, '%s',['#SBATCH --time=' num2str(Time) ':00:00' newline]);
        fprintf(fidout, '%s',['#SBATCH --no-requeue' newline]);
        fprintf(fidout, '%s',['#Setup of execution environment' newline]);
        fprintf(fidout, '%s',['#SBATCH --export=NONE' newline]);
        fprintf(fidout, '%s',['#SBATCH --get-user-env' newline]);
        fprintf(fidout, '%s',['#SBATCH --account=e05-surfin-clo' newline]);
        fprintf(fidout, '%s',['#SBATCH --partition=standard' newline]);
        fprintf(fidout, '%s',['#SBATCH --qos=lowpriority' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['module load slurm_setup' newline]);
        fprintf(fidout, '%s',['module load list' newline]);
        fprintf(fidout, '%s',['module load cp2k' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['# Ensure OMP_NUM_THREADS is consistent with cpus-per-task above' newline]);
        fprintf(fidout, '%s',['export OMP_NUM_THREADS=1' newline]);
        fprintf(fidout, '%s',['export OMP_PLACES=cores' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['cd ./' newline]);
        fprintf(fidout, '%s',['' newline]);
        fprintf(fidout, '%s',['srun --hint=nomultithread --distribution=block:block cp2k.psmp  ' OutFldr '-1.restart 2>&1 > out.log' newline]);
        fprintf(fidout, '%s',['' newline]);
        
        fclose(fidout);
    end
end