from aiida.engine import ToContext, WorkChain, calcfunction
from aiida.orm import AbstractCode, Int, Str
from aiida.plugins.factories import CalculationFactory, WorkflowFactory
from aiida.orm import load_computer
from aiida_shell import launch_shell_job
from aiida.engine import submit, run
from aiida.orm import SinglefileData
from io import StringIO
import os
import socket
import numpy as np
from spectra.spectra_vg_corr2 import spectrum_vg
from spectra.felixscript import make_output_file
from pathlib import Path

def engrad2np(natoms):
    grads = np.empty((4, natoms, 3))

    with open('final_step/grads.out', 'r') as f_in:
        for i, line in enumerate(f_in):
            if 3 <= i < (natoms + 3):
                grads[1][i-3] = [float(x) for x in line.split()[-3:]]
            if (natoms + 7) <= i < (natoms * 2 + 7):
                grads[2][i - natoms -7] = [float(x) for x in line.split()[-3:]]
            if (natoms * 2 + 11) <= i < (natoms * 3 + 11):
                grads[3][i - natoms * 2 - 11] = [float(x) for x in line.split()[-3:]]
            if (natoms * 3 + 15) <= i < (natoms * 4 + 15):
                grads[0][i - natoms * 3 - 15] = [float(x) for x in line.split()[-3:]]

    np.save('final_step/grads.npy', grads)

def mu2np():

    dips = np.zeros((6,6,3))
    excs = np.zeros(6)

    with open('final_step/ES_results.out','r') as f_in:
        for i, line in enumerate(f_in):
            if i in range(50, 55):
                dips[0][i-49] = [float(x) for x in line.split()[2:5]]
                dips[i-49][0] = [float(x) for x in line.split()[2:5]]
                excs[i-49] = float(line.split()[5])
            if i in range(59, 63):
                dips[1][i-57] = [float(x) for x in line.split()[2:5]]
                dips[i-57][1] = [float(x) for x in line.split()[2:5]]
            if i in range(64, 67):
                dips[2][i-61] = [float(x) for x in line.split()[2:5]]
                dips[i-61][2] = [float(x) for x in line.split()[2:5]]
            if i in range(68, 70):
                dips[3][i-64] = [float(x) for x in line.split()[2:5]]
                dips[i-64][3] = [float(x) for x in line.split()[2:5]]
            if i in range(71, 72):
                dips[4][i-66] = [float(x) for x in line.split()[2:5]]
                dips[i-66][4] = [float(x) for x in line.split()[2:5]]

        np.save('final_step/dipmoments.npy',dips[:4,:4,:3])
        np.save('final_step/excenergies.npy', excs[:4])


def communicate(message):
    # Prepare username and header and send them
    # We need to encode username to bytes, then count number of bytes and prepare header of fixed size, that we encode to bytes as well
    # message = (message).encode('utf-8')
    # message_header = f"{len(message):<{HEADER_LENGTH}}".encode('utf-8')
    # self.silasocket.send(message_header+ message)
    print(message)

class SilaLaserWorkChain(WorkChain):

    @classmethod
    def define(cls,spec):
        """Specify inputs and outputs."""
        super().define(spec)
        spec.input('smiles', valid_type=Str)
        spec.outline(
            cls.make_3d_struct,
            cls.laser_xtb_crest,
            cls.laser_orca_freq,
            cls.laser_orca_sp_nacsoc,
            cls.laser_orca_opt,
            cls.laser_orca_comb,
            cls.final_step,
        )
    def make_3d_struct(self):
        results,node = launch_shell_job(
            'obabel',
            arguments=[f"-:{self.inputs.smiles.value}", "--gen3d", "-o", "xyz"],
        )
        filecontent3d = results['stdout'].get_content()
        self.ctx.struct_3d = results['stdout'].get_content()

        communicate("openbabel success")
        return

    def laser_xtb_crest(self):
        communicate("beginning xtb crest job")
        results,node = launch_shell_job(
            '/home/a/aspuru/malcolms/laser/laser_xtb_crest.sh',
            arguments=["/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh", "{structfile}"],
            nodes={
                'structfile': SinglefileData(StringIO(self.ctx.struct_3d), filename="structfile.xyz")
            },
            outputs=[
                "crest", 
                "S0_XTB",
                "T1_XTB"
            ],
            metadata={
                'options': {
                    'computer': load_computer('niagara'),
                    'max_memory_kb':False,
                    'max_wallclock_seconds':10800,
                    'import_sys_environment': False,
                    'resources':{
                        'num_mpiprocs_per_machine':40,
                        'num_machines': 1
                    }                        
                }
            },
        )
        communicate("getting crest info..")
        # get crest information
        self.ctx.crest_conformers_xyz = results['crest'].get_object_content("crest_conformers.xyz")
        self.ctx.crest_best = results['crest'].get_object_content("crest_best.xyz")

        communicate("getting S0_XTB results...")
        # get S0 xtb results
        self.ctx.S0_XTB_xyz = results['S0_XTB'].get_object_content("xtbopt.xyz")
        # self.ctx.S0_XTB_out = results['S0_XTB'].get_object_content("98.out")
        self.ctx.S0_XTB_hessian = results['S0_XTB'].get_object_content("hessian")

        with open("xtb_results/S0_hessian", "w") as f:
            f.write(self.ctx.S0_XTB_hessian)
        with open("xtb_results/S0_XTB_xyz.xyz", "w") as f:
            f.write(self.ctx.S0_XTB_xyz)

        communicate("getting T1_XTB results...")
        # get T1 xtb results
        self.ctx.T1_XTB_xyz = results['T1_XTB'].get_object_content("xtbopt.xyz")
        # self.ctx.T1_XTB_out = results['T1_XTB'].get_object_content("98.out")
        self.ctx.T1_XTB_hessian = results['T1_XTB'].get_object_content("hessian")

        with open("xtb_results/T1_XTB_xyz.xyz", "w") as f:
            f.write(self.ctx.T1_XTB_xyz)
        with open("xtb_results/T1_hessian", "w") as f:
            f.write(self.ctx.T1_XTB_hessian)

        return
    
    def laser_orca_freq(self):
        communicate("starting orca_freq")
        results,node = launch_shell_job(
            '/home/a/aspuru/malcolms/laser/laser_orca_freq.sh',
            arguments=["/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh", "{structfile}"],
            nodes={
                'structfile': SinglefileData(StringIO(self.ctx.S0_XTB_xyz), filename="S0_XTB.xyz")
            },
            metadata={
                'options': {
                    'computer': load_computer('niagara'),
                    'import_sys_environment': False,
                    'max_wallclock_seconds': 10800, 
                    'resources':{
                        'num_mpiprocs_per_machine':40,
                        'num_machines': 4
                    }
                }
            },
            outputs=["freq"]
        )

        # get orca frequency job information
        self.ctx.orca_freq_xyz= results['freq'].get_object_content("FREQ_OPT.xyz")
        # self.ctx.orca_freq_hess= results['freq'].get_object_content("FREQ_OPT.hess")
        self.ctx.orca_freq_out= results['freq'].get_object_content("FREQ_OPT.out")

        with open("orca_freq_results/FREQ_OPT.xyz", "w") as f:
            f.write(self.ctx.orca_freq_xyz)
        # with open("orca_freq_results/FREQ_OPT.hess", "w") as f:
        #     f.write(self.ctx.orca_freq_hess)
        with open("orca_freq_results/FREQ_OPT.out", "w") as f:
            f.write(self.ctx.orca_freq_out)
        
        
        return
    
    def laser_orca_sp_nacsoc(self):
        communicate("starting orca_sp_nacsoc")
        results,node = launch_shell_job(
            '/home/a/aspuru/malcolms/laser/laser_orca_sp_nacsoc.sh',
            arguments=["/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh", "{structfile}"],
            nodes={
                'structfile': SinglefileData(StringIO(self.ctx.orca_freq_xyz), filename="orca_freq.xyz")
            },
            metadata={
                'options': {
                    'computer': load_computer('niagara'),
                    'import_sys_environment': False,
                    'max_wallclock_seconds': 10800,
                    'resources':{
                        'num_mpiprocs_per_machine':40,
                        'num_machines': 1
                    }
                }
            },
            outputs=["S0_tdaNAC", "S0_tdaSOC"]
        )


        # get orca_sp_nacsoc information
        self.ctx.S0_tdaNAC_out= results['S0_tdaNAC'].get_object_content("S0_TDA_SP.out")
        self.ctx.S0_tdaSOC_out= results['S0_tdaSOC'].get_object_content("S0_TDA_SP.out")
        self.ctx.S0_tdaSOC_inp= results['S0_tdaSOC'].get_object_content("S0_TDA_SP.inp")
        self.ctx.S0_tdaNAC_inp= results['S0_tdaNAC'].get_object_content("S0_TDA_SP.inp")
        self.ctx.S0_tdaNAC_molden= results['S0_tdaNAC'].get_object_content("S0_TDA_SP.molden")
        self.ctx.S0_tdaSOC_molden= results['S0_tdaSOC'].get_object_content("S0_TDA_SP.molden")


        with open("orca_sp_nacsoc_results/S0_tdaNAC_SP.out", "w") as f:
            f.write(self.ctx.S0_tdaNAC_out)
        with open("orca_sp_nacsoc_results/S0_tdaSOC_SP.out", "w") as f:
            f.write(self.ctx.S0_tdaSOC_out)
        with open("orca_sp_nacsoc_results/S0_tdaSOC_SP.inp", "w") as f:
            f.write(self.ctx.S0_tdaSOC_inp)
        with open("orca_sp_nacsoc_results/S0_tdaNAC_SP.molden", "w") as f:
            f.write(self.ctx.S0_tdaNAC_molden)
        with open("orca_sp_nacsoc_results/S0_tdaSOC_SP.molden", "w") as f:
            f.write(self.ctx.S0_tdaSOC_molden)
        

    def laser_orca_opt(self):
        communicate("starting orca opt job")
        results,node = launch_shell_job(
            '/home/a/aspuru/malcolms/laser/run_orca_opts.sh',
            arguments=["/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh"],
            nodes={
                'S0_XTB_OPT_XYZ':SinglefileData(StringIO(self.ctx.S0_XTB_xyz), filename="S0_XTB_OPT_XYZ.xyz"),
                'T1_XTB_OPT_XYZ':SinglefileData(StringIO(self.ctx.T1_XTB_xyz), filename="T1_XTB_OPT_XYZ.xyz"),
            },
            metadata={
                'options': {
                    'computer': load_computer('niagara'),
                    'import_sys_environment': False,
                    'max_wallclock_seconds': 28800,
                    'resources':{
                        'num_mpiprocs_per_machine':40,
                        'num_machines': 1
                    }
                }
            },
            outputs = ["S1", "T1", "T2", "T3", "T4"]
        )

        
        self.ctx.S1 = results['S1'].get_object_content("S1_OPT.xyz")
        self.ctx.T1 = results['T1'].get_object_content("T1_OPT.xyz")
        self.ctx.T2 = results['T2'].get_object_content("T2_OPT.xyz")
        self.ctx.T3 = results['T3'].get_object_content("T3_OPT.xyz")
        self.ctx.T4 = results['T4'].get_object_content("T4_OPT.xyz")

        with open("orca_opt_results/S1_OPT.xyz", "w") as f:
            f.write(self.ctx.S1)
        with open("orca_opt_results/T1_OPT.xyz", "w") as f:
            f.write(self.ctx.T1)
        with open("orca_opt_results/T2_OPT.xyz", "w") as f:
            f.write(self.ctx.T2)
        with open("orca_opt_results/T3_OPT.xyz", "w") as f:
            f.write(self.ctx.T3)
        with open("orca_opt_results/T4_OPT.xyz", "w") as f:
            f.write(self.ctx.T4)
        return

    def laser_orca_comb(self):
        communicate(" beginning comb step")
        outputs =[
            "S1_OPT_Singlet",  
            "S1_OPT_Triplet",
            "T1_OPT_Singlet",  
            "T1_OPT_Triplet",  
            "T2_OPT_Singlet",
            "T2_OPT_Triplet", 
            "T3_OPT_Singlet",
            "T3_OPT_Triplet",
            "T4_OPT_Singlet",  
            "T4_OPT_Triplet"
            ]
        results,node = launch_shell_job(
            '/home/a/aspuru/malcolms/laser/laser_orca_comb.sh',
            arguments=["/home/a/aspuru/malcolms/laser/testing/lib/funcs.sh", "{S1}", "{T1}", "{T2}", "{T3}", "{T4}" ],
            nodes={
                'S1':SinglefileData(StringIO(self.ctx.S1), filename="S1_OPT.xyz"),
                'T1':SinglefileData(StringIO(self.ctx.T1), filename="T1_OPT.xyz"),
                'T2':SinglefileData(StringIO(self.ctx.T2), filename="T2_OPT.xyz"),
                'T3':SinglefileData(StringIO(self.ctx.T3), filename="T3_OPT.xyz"),
                'T4':SinglefileData(StringIO(self.ctx.T4), filename="T4_OPT.xyz"),
            },
            metadata={
                'options': {
                    'computer': load_computer('niagara'),
                    'import_sys_environment': False,
                    'max_wallclock_seconds': 28800,
                    # 'queue_name': "debug", 
                    'resources':{
                        'num_mpiprocs_per_machine':40,
                        'num_machines': 1
                    }
                }
            },
            outputs=outputs,
        )

        for dir in outputs:

            with open(f"orca_comb_results/{dir}/TDA_SP_COMB.out", "w" ) as f:
                f.write(results[dir].get_object_content("TDA_SP_COMB.out"))

            with open(f"orca_comb_results/{dir}/TDA_SP_COMB.inp", "w" ) as f:
                f.write(results[dir].get_object_content("TDA_SP_COMB.inp"))

            with open(f"orca_comb_results/{dir}/TDA_SP_COMB_property.txt", "w" ) as f:
                f.write(results[dir].get_object_content("TDA_SP_COMB_property.txt"))
        
        
        return
    
    def final_step(self):

        print("Executing final step")

        with open("xtb_results/S0_XTB_xyz.xyz", "r") as f:
            natoms = int(f.readline())


        hess_arr = np.empty(pow(natoms*3, 2))

        with open('xtb_results/S0_hessian', 'r') as f_in:
            count = 0
            for i, line in enumerate(f_in):
                if i == 0:
                    continue
                vibs = [float(x) for x in line.split()]
                for vib in vibs:
                    hess_arr[count] = vib
                    count += 1
        np.save('final_step/hessian.npy', hess_arr.reshape((natoms*3, natoms*3)))

    
        results,node = launch_shell_job(
            "/home/malcolm/sila-aiida/final_step_aiida.sh",
            nodes={
                'NACSPout': SinglefileData(StringIO(Path('./orca_sp_nacsoc_results/S0_tdaNAC_SP.out').read_text()), filename="S0_tdaNAC_SP.out"),
                'S0XTBxyz': SinglefileData(StringIO(Path('./xtb_results/S0_XTB_xyz.xyz').read_text()), filename="S0_XTB_xyz.xyz"),
                'molden': SinglefileData(StringIO(Path('./orca_sp_nacsoc_results/S0_tdaNAC_SP.molden').read_text()), filename="S0_tdaNAC_SP.molden"),
                'ESinfo': SinglefileData(StringIO(Path('./spectra/ES_info.txt').read_text()), filename="ES_info.txt"),


            },
            arguments=[str(natoms)],
            outputs=[
                "ES_results.out", 
                "grads.out",
                "orca_mo.out",
                "ES_info.out"
            ],
        )

        print(results.keys())

        print(results["stderr"].get_content())

        with open("final_step/ES_results.out", "w") as f:
            f.write(results["ES_results_out"].get_content())
        with open("final_step/grads.out", "w") as f:
            f.write(results["grads_out"].get_content())
        with open("final_step/orca_mo.out", "w") as f:
            f.write(results["orca_mo_out"].get_content())
        
        with open("final_step/ES_info.out", "w") as f:
            f.write(results["ES_info_out"].get_content())
        
        engrad2np(natoms)

        mu2np()

        spectrum_vg()
        
        make_output_file()
    


        





   



        





   
