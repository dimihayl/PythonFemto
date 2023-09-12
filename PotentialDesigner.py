import os,sys
import re
import ROOT
import math
import optuna
import subprocess
import time
import shutil

CPU_HOURS = 24

script_name = "/home/dimihayl/Software/LocalFemto/bin/PotentialDesignerEngine"

#sinfo: kta,long,xtralong
SBATCH_BASE = 'sbatch --partition=long'
SBATCH_BASE += ' --mem-per-cpu=1000'
SBATCH_BASE += ' --time='+str(CPU_HOURS*60)
SBATCH_BASE += ' --job-name=PotDes'
#SBATCH_BASE += ' --constraint=INTEL'
SBATCH_BASE += ' /home/ktas/gu47det/DimiSoftware/UsmaniLayeredFitAnalysis/UsmaniLayeredFit_Fitter.sh'

FmToNu = 5.067731237e-3

#studyR = optuna.create_study(direction='minimize', sampler=optuna.samplers.RandomSampler())
studyT = optuna.create_study(direction='minimize', sampler=optuna.samplers.TPESampler())

par1_range = [0, 0]
par2_range = [0, 0]
par3_range = [0, 0]
par4_range = [0, 0]

#par[0] are the number of pars in the ERE
def fit_scattering_pars(x, par):
    MOM = x[0]
    Npars = round(par[0])
    #print(Npars)
    if Npars==0:
        return 0
    ERE = 0
    if Npars>=1:
        ERE += 1./(par[1]*FmToNu)
    if Npars>=2:
        ERE += 0.5*MOM*MOM*par[2]*FmToNu
    if Npars>=3:
        for iPar in range(2,Npars):
            ERE += par[iPar+1] * (FmToNu**(iPar*2-1)) * (MOM**(iPar*2))
    return math.atan(MOM/ERE)

# Target function to be optimized by Optuna
def target_function(trial):
    # Sample the parameters within the specified ranges
    v_par1 = trial.suggest_float('v_par1', *par1_range)
    v_par2 = trial.suggest_float('v_par2', *par2_range)
    v_par3 = trial.suggest_float('v_par3', *par3_range)
    v_par4 = trial.suggest_float('v_par4', *par4_range)
    return 0

def run_script(args):
    subprocess.Popen([script_name] + [args], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def parse_float_range(range_str):
    return [float(value) for value in range_str.split()]

def get_yes_no_input(prompt):
    while True:
        user_input = input(prompt).strip().lower()
        if user_input in ["yes", "y"]:
            return True
        elif user_input in ["no", "n"]:
            return False
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")

def main():
    # Check if the correct number of arguments are provided
    if len(sys.argv) != 2:
        print("Usage: python script.py <filename>")
        return

    BaseFileName = sys.argv[1]
    filename_arg = BaseFileName+'.txt'
    #ROOT_FILE_NAME = sys.argv[1]+'.root'

    f_goal = 0
    f_err = 0
    d_goal = 0
    d_err = 0
    #par1 = 0
    #par2 = 0
    #par3 = 0
    #par4 = 0
    pot = ""
    M1 = 0
    M2 = 0
    PW = 0
    q1q2 = 0
    #a root file containing phase shifts to be fitted
    PhaseShifts = ""
    kMin = 0
    kMax = 100
    kBin = 50
    CPU = 1
    FIT_TIMEOUT = 1#in minutes
    TIMEOUT = 20#in minutes
    WakeUp = 2.5#in seconds
    nPars = 3

    if not os.path.exists(filename_arg):
        sys.exit(f"The file '{filename_arg}' does not exist.")

    # Open the file for reading
    with open(filename_arg, 'r') as file:
        for line in file:
            single_line = line.strip()
            line_split = single_line.split()
            #print(line_split)
            if line_split[0].lower()=='f_goal':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of f_goal')
                try:
                    f_goal = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: f_goal must be a number.")
            if line_split[0].lower()=='f_err':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of f_err')
                try:
                    f_err = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: f_err must be a number.")
            if line_split[0].lower()=='d_goal':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of d_goal')
                try:
                    d_goal = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: d_goal must be a number.")
            if line_split[0].lower()=='d_err':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of d_err')
                try:
                    d_err = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: d_err must be a number.")
            if line_split[0].lower()=='par1':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par1')
                try:
                    par1_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par1 must be a number.")
                if len(line_split)==3:
                    try:
                        par1_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par1 must be a number.")
                else:
                    par1_range[1] = par1_range[0]
            if line_split[0].lower()=='par2':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par2')
                try:
                    par2_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par2 must be a number.")
                if len(line_split)==3:
                    try:
                        par2_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par2 must be a number.")
                else:
                    par2_range[1] = par2_range[0]
            if line_split[0].lower()=='par3':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par3')
                try:
                    par3_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par3 must be a number.")
                if len(line_split)==3:
                    try:
                        par3_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par3 must be a number.")
                else:
                    par3_range[1] = par3_range[0]
            if line_split[0].lower()=='par4':
                if len(line_split)!=2 and len(line_split)!=3:
                    sys.exit('Wrong set up of par4')
                try:
                    par4_range[0] = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: par4 must be a number.")
                if len(line_split)==3:
                    try:
                        par4_range[1] = float(line_split[2])  # Try converting the string to a float
                    except ValueError:
                        sys.exit("Error: par4 must be a number.")
                else:
                    par4_range[1] = par4_range[0]
            if line_split[0].lower()=='m1':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of M1')
                try:
                    M1 = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: M1 must be a number.")
                if M1<=0:
                    sys.exit("Error: M1 must be non-zero and positive.")
            if line_split[0].lower()=='m2':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of M2')
                try:
                    M2 = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: M2 must be a number.")
                if M2<=0:
                    sys.exit("Error: M2 must be non-zero and positive.")
            if line_split[0].lower()=='kmin':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kMin')
                try:
                    kMin = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: kMin must be a number.")
                if kMin<0:
                    sys.exit("Error: kMin must be non-negative.")
            if line_split[0].lower()=='kmax':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kMax')
                try:
                    kMax = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: kMax must be a number.")
                if kMax<=0:
                    sys.exit("Error: kMax must be positive.")
            if line_split[0].lower()=='kbin':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of kMax')
                if not line_split[1].isdigit():
                    sys.exit('kBin is not a digit')
                kBin = int(line_split[1])
                if kBin<=0:
                    sys.exit("Error: kBin must be positive.")
            if line_split[0].lower()=='eps':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of eps')
                try:
                    eps = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: eps must be a number.")
                if eps<1e-10 or eps>1e-5:
                    sys.exit("Error: eps must be >1e-10 and <1e-5 for CATS to properly work.")
            if line_split[0].lower()=='pot':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of pot')
                pot = line_split[1]
                if pot not in ["Gauss","DoubleGauss","Yukawa","YukawaDLM"]:
                    sys.exit("Error: pot should be selected from 'Gauss, DoubleGauss, Yukawa, YukawaDLM'")
            if line_split[0].lower()=='pw':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of pw')
                pw_str = line_split[1]
                if pw_str.lower() not in ["s","p","d","f"]:
                    sys.exit("Error: pw should be selected from 's, p, d, f'")
                if pw_str.lower()=='s':
                    PW = 0
                if pw_str.lower()=='p':
                    PW = 1
                if pw_str.lower()=='d':
                    PW = 2
                if pw_str.lower()=='f':
                    PW = 3
            if line_split[0].lower()=='coulomb':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of Coulomb')
                if not line_split[1].isdigit():
                    sys.exit('Coulomb is not a digit')
                q1q2 = int(line_split[1])
            if line_split[0].lower()=='phaseshifts':
                if len(line_split)!=3:
                    sys.exit('Wrong set up of PhaseShifts')
                PhaseShifts = line_split[1]
                PhaseShiftsHisto = line_split[2]
            if line_split[0].lower()=='cpu':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of CPU')
                if not line_split[1].isdigit():
                    sys.exit('CPU is not a digit')
                CPU = int(line_split[1])
                if CPU<=0:
                    sys.exit("Error: CPU must be positive.")
                if CPU>32:
                    sys.exit("Error: CPU must be <= 32.")
            if line_split[0].lower()=='npars':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of nPars')
                if not line_split[1].isdigit():
                    sys.exit('nPars is not a digit')
                nPars = int(line_split[1])
                if nPars<=0:
                    sys.exit("Error: nPars must be positive.")
                if nPars>8:
                    sys.exit("Error: nPars must be <= 8.")

            if line_split[0].lower()=='timeout':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of TIMEOUT')
                if not line_split[1].isdigit():
                    sys.exit('TIMEOUT is not a digit')
                TIMEOUT = int(line_split[1])
                if TIMEOUT<=0:
                    sys.exit("Error: TIMEOUT must be positive.")
            if line_split[0].lower()=='fit_timeout':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of FIT_TIMEOUT')
                if not line_split[1].isdigit():
                    sys.exit('FIT_TIMEOUT is not a digit')
                FIT_TIMEOUT = int(line_split[1])
                if FIT_TIMEOUT<=0:
                    sys.exit("Error: FIT_TIMEOUT must be positive.")
            if line_split[0].lower()=='wakeup':
                if len(line_split)!=2:
                    sys.exit('Wrong set up of WakeUp')
                try:
                    WakeUp = float(line_split[1])  # Try converting the string to a float
                except ValueError:
                    sys.exit("Error: WakeUp must be a number.")

        if par1_range[0]>par1_range[1]:
            sys.exit('Error: par1 minimum value is larger than the maximum.')
        if par2_range[0]>par2_range[1]:
            sys.exit('Error: par2 minimum value is larger than the maximum.')
        if par3_range[0]>par3_range[1]:
            sys.exit('Error: par3 minimum value is larger than the maximum.')
        if par4_range[0]>par4_range[1]:
            sys.exit('Error: par4 minimum value is larger than the maximum.')

        if par2_range[0] <= 0:
            sys.exit('par2 is the width of the potential and should be non-zero and positive.')
        if par4_range[0] <= 0 and pot=="DoubleGauss" :
            sys.exit('par4 is the width of the potential and should be non-zero and positive.')

        if pot not in ["Gauss","DoubleGauss","Yukawa","YukawaDLM"]:
            sys.exit("Error: pot should be selected from 'Gauss, DoubleGauss, Yukawa, YukawaDLM'")
        if M1<=0 or M2<=0:
            sys.exit("Error: Both M1 and M2 should be defined and be strictly positive.")
        if f_goal!=0 and PhaseShifts != "":
            sys.exit("Error: The input should be defined either through f_goal or a PhaseShifts, but not both!")

        if kMin<0:
            sys.exit("Error: kMin cannot be negative!")
        if kMin>100:
            sys.exit("Error: kMin should be small, certainly not above 100 MeV!")
        if kMin>=kMax:
            sys.exit("Error: kMin should be smaller than kMax!")
        if kMax<0:
            sys.exit("Error: kMax cannot be negative!")
        if kBin<=0:
            sys.exit("Error: kBin should be positive!")

        if f_err==0:
            f_err = abs(0.01*f_goal)
            if f_err < 0.01:
                f_err = 0.01
        if d_err==0:
            d_err = abs(0.01*d_goal)
            if d_err < 0.01:
                d_err = 0.01



        if PhaseShifts!="":
            if not os.path.exists(PhaseShifts):
                sys.exit(f"The file '{PhaseShifts}' does not exist.")
            # Open a ROOT file
            ps_file = ROOT.TFile.Open(PhaseShifts)

            # Check if the file was opened successfully
            if not ps_file or ps_file.IsZombie():
                sys.exit(f"The file '{PhaseShifts}' cannot be opened.")

            PS_TO_FIT = ps_file.Get(PhaseShiftsHisto)
            if not PS_TO_FIT:
                sys.exit(f"The histogram: '{PhaseShiftsHisto}' could not be retrieved")

            # Close the ROOT file when you're done
            ps_file.Close()

        elif PW==0:
            PhaseShifts = 'FIT_F0_D0_DIRECTLY'
            PS_TO_FIT = ROOT.TH1F("PS_TO_FIT","PS_TO_FIT",kBin,kMin,kMax)
            for uMom in range(0,kBin):
                MOM = PS_TO_FIT.GetBinCenter(uMom+1)
                ps_val = math.atan(MOM/(1./(f_goal*FmToNu)+0.5*MOM*MOM*d_goal*FmToNu))
                PS_TO_FIT.SetBinContent(uMom+1,ps_val)
        else:
            sys.exit("ERROR: To study l!=0 on has to provide the phase shifts as a direct input!")
        PS_TO_FIT.SetName("PS_TO_FIT")

        Graph_PS_TO_FIT = ROOT.TGraph()
        Graph_PS_TO_FIT.SetName("Graph_PS_TO_FIT")
        for uMom in range(0,kBin):
            MOM = PS_TO_FIT.GetBinCenter(uMom+1)
            Graph_PS_TO_FIT.SetPoint(uMom,MOM,PS_TO_FIT.GetBinContent(uMom+1))

        #ROOT_FILE = ROOT.TFile.Open(ROOT_FILE_NAME,'recreate')
        #PS_TO_FIT.Write()

        f_err_current = f_err*10
        d_err_current = d_err*10
        UniqueID = 0
        BestEstimator = 1e64
        Best_f = 1e64
        Best_d = 1e64
        Best_par1 = 0
        Best_par2 = 0
        Best_par3 = 0
        Best_par4 = 0
        Progress_f = 0
        Progress_d = 0

        TotExeTime = 0
        StartTime = time.time()/60.
        #we re-iterate until we reach our convergence criteria
        while (f_err_current>f_err or d_err_current>d_err) and TotExeTime<TIMEOUT:
            TrialsCPU = []
            FilesIn = []
            FilesOut = []
            CurrentPar1 = []
            CurrentPar2 = []
            CurrentPar3 = []
            CurrentPar4 = []
            for iCPU in range(0,CPU):
                UniqueID += 1
                fileBaseJob = BaseFileName+"_"+str(UniqueID)
                fileNameIn = fileBaseJob+".fit"
                fileNameOut = fileBaseJob+".root"
                FilesIn.append(fileNameIn)
                FilesOut.append(fileNameOut)
                #EstimatorsCPU.append(0)

                trial = studyT.ask()
                TrialsCPU.append(trial)
                target_function(trial)


                try:
                    shutil.copy(filename_arg, fileNameIn)
                except FileNotFoundError:
                    print(f"Source file '{filename_arg}' not found.")
                except Exception as e:
                    sys.exit(f"An error occurred: {e}")


                # Check if the file exists
                if os.path.exists(fileNameOut):
                    try:
                        # Delete the file
                        os.remove(fileNameOut)
                    except OSError:
                        pass  # Ignore errors silently if they occur during file deletion

                with open(fileNameIn, 'a') as fjob:
                    fjob.write( "UniqueID   "+str(UniqueID)+"\n")
                    fjob.write( "v_par1     "+str(trial.params['v_par1'])+"\n")
                    fjob.write( "v_par2     "+str(trial.params['v_par2'])+"\n")
                    fjob.write( "v_par3     "+str(trial.params['v_par3'])+"\n")
                    fjob.write( "v_par4     "+str(trial.params['v_par4'])+"\n")
                    CurrentPar1.append(trial.params['v_par1'])
                    CurrentPar2.append(trial.params['v_par2'])
                    CurrentPar3.append(trial.params['v_par3'])
                    CurrentPar4.append(trial.params['v_par4'])

                run_script(fileBaseJob)
            #for iCPU in range(0,CPU)

            #wait for all jobs to finish
            StartWaitTime = time.time()/60.
            NumRunningJobs = CPU
            while NumRunningJobs>0 and (time.time()/60.-StartWaitTime)<FIT_TIMEOUT:
                NumRunningJobs = CPU
                for iCPU in range(0,CPU):
                    if os.path.exists(FilesOut[iCPU]):
                        try:
                            NumRunningJobs -= 1
                            time.sleep(WakeUp*0.1)
                        except OSError:
                            time.sleep(WakeUp)
            time.sleep(0.5)
            if NumRunningJobs>0:
                print('WARNING: Some jobs have failed')
            #after all jobs are done, we collect the output
            for iCPU in range(0,CPU):
                if os.path.exists(FilesOut[iCPU]):
                    try:
                        root_file = ROOT.TFile.Open(FilesOut[iCPU])
                        hPhaseShifts = root_file.Get("hPhaseShifts")
                        hPotential = root_file.Get("hPotential")
                        Estimator = 0.0
                        if PhaseShifts=='FIT_F0_D0_DIRECTLY':
                            fit_ps_def = ROOT.TF1("fit_ps_def",fit_scattering_pars,kMin,kMax,nPars+1)
                            fit_ps_def.FixParameter(0,nPars)
                            fit_ps_def.SetParameter(1,f_goal)
                            fit_ps_def.SetParameter(2,d_goal)
                            for iPar in range(2,nPars):
                                fit_ps_def.SetParameter(iPar+1,0)
                            for uMom in range(0,kBin):
                                MOM = hPhaseShifts.GetBinCenter(uMom+1)
                                hPhaseShifts.SetBinError(uMom+1,0.001)
                            hPhaseShifts.Fit(fit_ps_def, "Q, S, N, R, M")
                            Current_f = fit_ps_def.GetParameter(1)
                            Current_d = fit_ps_def.GetParameter(2)
                            Estimator += ( (Current_f-f_goal)/f_err )**2
                            Estimator += ( (Current_d-d_goal)/d_err )**2
                        else:
                            for uMom in range(0,kBin):
                                MOM = hPhaseShifts.GetBinCenter(uMom+1)
                                #a dummy error that has linear scaling. This is to consider that we have extra terms to ERE at large k*
                                ERR = 0.01*MOM
                                kCotGoal = MOM/math.tan(hPhaseShifts.GetBinContent(uMom+1))
                                kCotFit = MOM/math.tan(Graph_PS_TO_FIT.Eval(MOM))
                                #Estimator += (hPhaseShifts.GetBinContent(uMom+1)-Graph_PS_TO_FIT.Eval(MOM))**2
                                Estimator += ( (kCotGoal-kCotFit)/ERR )**2
                        if BestEstimator>Estimator:
                            #BestTrial = TrialsCPU[iCPU]
                            for uMom in range(0,kBin):
                                hPhaseShifts.SetBinError(uMom+1,0.001)
                            fit_ps = ROOT.TF1("fit_ps",fit_scattering_pars,kMin,kMax,nPars+1)
                            fit_ps.FixParameter(0,nPars)
                            fit_ps.SetParameter(1,f_goal)
                            fit_ps.SetParameter(2,d_goal)
                            for iPar in range(2,nPars):
                                fit_ps.SetParameter(iPar+1,0)
                            #fit_ps = ROOT.TF1("fit_ps","TMath::ATan(x/(1./[0]/5.067731237e-3+0.5*x*x*[1]*5.067731237e-3))",kMin,kMax)
                            #fit_ps.SetParameter(0,f_goal)
                            #fit_ps.SetParameter(1,d_goal)
                            hPhaseShifts.Fit(fit_ps, "Q, S, N, R, M")
                            Best_f = fit_ps.GetParameter(1)
                            Best_d = fit_ps.GetParameter(2)

                            BestEstimator = Estimator
                            Best_par1 = CurrentPar1[iCPU]
                            Best_par2 = CurrentPar2[iCPU]
                            Best_par3 = CurrentPar3[iCPU]
                            Best_par4 = CurrentPar4[iCPU]
                            f_err_current = abs(Best_f-f_goal)
                            d_err_current = abs(Best_d-d_goal)
                            Progress_f = f_err/f_err_current
                            Progress_d = d_err/d_err_current
                            # The message you want to update
                            if Progress_f>1:
                                Progress_f = 1
                            if Progress_d>1:
                                Progress_d = 1
                            print(f"Progress: f({Progress_f*100:.1f} %) = {Best_f:.3f} fm     d({Progress_d*100:.1f} %) = {Best_d:.3f} fm        ", end="\r")
                            best_file = ROOT.TFile.Open(BaseFileName+"_BestSolution.root","recreate")
                            Graph_PS_TO_FIT.Write()
                            hPhaseShifts.Write()
                            hPotential.Write()
                            root_file.cd()
                            shutil.copy(FilesIn[iCPU], BaseFileName+"_BestSolution.txt")
                            with open(BaseFileName+"_BestSolution.txt", 'a') as fupdate:
                                fupdate.write( "f_best     "+str(Best_f)+"\n")
                                fupdate.write( "d_best     "+str(Best_d)+"\n")

                        studyT.tell(TrialsCPU[iCPU],Estimator)
                        root_file.Close()
                        #clean up
                        os.remove(FilesOut[iCPU])
                        os.remove(FilesIn[iCPU])
                    except OSError:
                        pass
            #for iCPU in range(0,CPU)
            TotExeTime = time.time()/60.-StartTime
        #while f_err_current>f_err or d_err_current>d_err

        print("\nSUMMARY:")
        if f_err_current<=f_err and d_err_current<=d_err:
            print(" - The desired convergence was reached!")
        else:
            print(" - The script terminated due to timeout!")
        print(f" - Obtained f = {Best_f:.3f} fm")
        print(f" - Obtained d = {Best_d:.3f} fm")
        print(" - The potential was "+pot)
        if pot=="DoubleGauss":
            print(f"   V1 = {Best_par1:.2f} MeV")
            print(f"   μ1 = {Best_par2:.5f} fm")
            print(f"   V2 = {Best_par3:.2f} MeV")
            print(f"   μ2 = {Best_par4:.5f} MeV")
        else:
            print(f"   V = {Best_par1:.2f} MeV")
            print(f"   μ = {Best_par2:.5f} fm")
        print(" - Detailed output in the "+BaseFileName+"_BestSolution"+".txt/root files.")


if __name__ == "__main__":
    main()
