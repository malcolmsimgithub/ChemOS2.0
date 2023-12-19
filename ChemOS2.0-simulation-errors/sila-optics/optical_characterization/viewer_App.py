import numpy as np
import matplotlib.pyplot as plt
import json
import pickle as pkl
import os
from tkinter import filedialog
import tkinter as tk
from tkinter import messagebox, ttk
import copy

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import time

from viewer_utils.spec_to_color.spec_to_color import Spec_to_Color
from viewer_utils.rader_chart import radar_factory
from viewer_utils.viewer_config import config_dict
from viewer_utils.heatmap import heatmap
# from classes import UV, Absorption, Emission


filepath = os.path.dirname(os.path.abspath(__file__))
image_file = '%s/temp/mol_image.png' %filepath


nF = 1.5 #refractive index of a film for gain calculation



class Application(tk.Frame):

    def __init__(self, master = None):
        super().__init__(master, height = 60, width = 230)
        self.pack()

        self.config = config_dict

        #self.master.geometry("300x300")
        self.master.title("Viewer")

        ### Menu
        menu_bar = tk.Menu(self, tearoff=0)
        # File
        menu_file = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=menu_file,  underline=0)
        menu_file.add_command(label="Open", command=self.file_open, underline=0, accelerator = 'Ctrl-O')
        menu_file.add_command(label="Exit", command=self.exit, underline=0, accelerator = 'Ctrl-Q')
        # File info
        # menu_info = tk.Menu(menu_bar, tearoff=0)
        # menu_bar.add_cascade(label="File_info", menu=menu_info,  underline=0)
        # # menu_info.add_command(label="Injection_info", command=self.show_injection_info, underline=0)
        # # menu_info.add_command(label="HPLC_info", command=self.show_HPLC_info, underline=0)
        # menu_info.add_command(label="Analysis_info", command=self.show_analysis_info, underline=0)
        self.master.config(menu=menu_bar)

        self.r2_val = tk.IntVar()
        self.text_filename = tk.StringVar()
        self.text_samplename = tk.StringVar()
        self.text_solvent = tk.StringVar()
        self.text_absmax = tk.StringVar()
        self.text_absend = tk.StringVar()
        self.text_PLmax = tk.StringVar()
        self.text_gainmax = tk.StringVar()
        self.text_relativeQY = tk.StringVar()
        self.text_PLlifetime = tk.StringVar()
        self.text_radiative_decay = tk.StringVar()
        self.text_PLretention = tk.StringVar()
        self.text_Absretention = tk.StringVar()
        self.text_uvref  = tk.StringVar()
        self.text_uvabs = tk.StringVar()
        self.text_gain_cs = tk.StringVar()
        self.text_gcs_PLretention = tk.StringVar()

        self.data = None
        self.r_buttons = None

        self.master.protocol("WM_DELETE_WINDOW", self.exit)

    #wrapper
    def data_method(func):
        def wrapper(self, *args, **kwargs):
            if not self.data:
                messagebox.showwarning('Warning', 'Please load results file')
            else: 
                func(self, *args, **kwargs)
        return wrapper   


    def open_graph_window(self):

        def on_closing():
            plt.close()
            self.win_graph.destroy()

        c = len(self.filenames)//36

        self.win_graph = tk.Toplevel()
        self.win_graph.wm_geometry('%sx900' %(1500 + 150 *c))
        self.win_graph.title('Optical Characterization results') 

        menu_bar = tk.Menu(self.win_graph, tearoff=0)

        menu_info = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File_info", menu=menu_info,  underline=0)
        menu_info.add_command(label="measurement_info", command=self.show_measure_info, underline=0)
        menu_info.add_command(label="job_info", command=self.show_job_info, underline=0)

        menu_plot = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="Plots", menu=menu_plot,  underline=0)
        menu_plot.add_command(label="Heat map", command=self.show_heatmap, underline=0)
        menu_plot.add_command(label="Scatter plot ", command= lambda : self.show_scatter_plots(False), underline=0)
        menu_plot.add_command(label="Scatter plot (with label)", command=lambda : self.show_scatter_plots(True), underline=0)

        self.win_graph.config(menu=menu_bar)

        self.create_widgets() 

        self.win_graph.protocol("WM_DELETE_WINDOW", on_closing)


    def create_widgets(self):

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        self.canvas0 = tk.Canvas(self.win_graph, bg = 'black', height = 280, width = 320)
        self.canvas0.place(x = 40, y = 40)

        ###Plot spectra
        self.fig1 = plt.figure(figsize=(7,7))
        self.ax1 = self.fig1.add_subplot(2,2,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance') #Absorption/PL 
        self.ax2 = self.ax1.twinx()
        self.ax2.set_ylabel('PL emission/energy nm-1 s-1')
        self.ax3 = self.fig1.add_subplot(2,2,2, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
        self.ax4 = self.ax3.twinx()
        self.ax4.set_ylabel('Gain factor/1E24 cm2 s')
        self.ax5 = self.fig1.add_subplot(2,2,3, xlabel = 'time / ns', ylabel = 'photon count')
        self.ax6 = self.fig1.add_subplot(2,2,4, xlabel = 'time / s', ylabel = 'Normalized Absorbance(365nm)/PL intensity')

        self.fig1.tight_layout()

        self.canvas1 = FigureCanvasTkAgg(self.fig1, master= self.win_graph)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().place(x=400, y = 30, height=650, width=950)

        ###Create chart
        self.spoke_labels = ['{}\n({:.1e}{})'.format(val['Label'], val['Normalize_at'], val['unit']) for val in self.config.values()]
        self.N = len(self.spoke_labels)
        self.theta = radar_factory(self.N, frame='polygon')      
        
        self.fig2 = plt.figure(figsize =(4,4))
        self.ax7 = self.fig2.add_subplot(1,1,1, projection = 'radar')

        self.ax7.set_rgrids([0.2, 0.4, 0.6, 0.8, 1.0], labels = ['','','','',''])
        self.ax7.set_varlabels(self.spoke_labels)

        self.fig2.tight_layout()

        self.canvas2 = FigureCanvasTkAgg(self.fig2, master= self.win_graph)
        self.canvas2.draw()
        self.canvas2.get_tk_widget().place(x=40, y = 350, height=330, width=340)
        


        x_offset = 40
        y_offset = 750
        #filename
        label_0 = tk.Label(self.win_graph, text = 'filename :', anchor = 'w').place(x=x_offset, y = y_offset - 30)
        label_1 = tk.Label(self.win_graph, text = 'sample name :', anchor = 'w').place(x = x_offset, y = y_offset)
        label_2 = tk.Label(self.win_graph, text = 'solvent :', anchor = 'w').place(x = x_offset, y = y_offset + 30)
        label_3 = tk.Label(self.win_graph, text = 'Absorption max :', anchor = 'w').place(x = x_offset, y = y_offset + 60)
        label_4 = tk.Label(self.win_graph, text = 'Absorption end :', anchor = 'w').place(x = x_offset, y = y_offset + 90)
        label_5 = tk.Label(self.win_graph, text = 'PL max :', anchor = 'w').place(x = x_offset + 350, y = y_offset - 30)
        label_6 = tk.Label(self.win_graph, text = 'Gain factor max(cm2 s) :', anchor = 'w').place(x = x_offset + 350, y = y_offset)
        label_7 = tk.Label(self.win_graph, text = 'Relative QY :', anchor = 'w').place(x = x_offset + 350, y = y_offset + 30)
        label_8 = tk.Label(self.win_graph, text = 'PL lifetime :', anchor = 'w').place(x = x_offset + 350, y = y_offset + 60)
        label_9 = tk.Label(self.win_graph, text = 'radiative decay constant:', anchor = 'w').place(x = x_offset + 350, y = y_offset + 90)
        label_10 = tk.Label(self.win_graph, text = 'PL retention :', anchor = 'w').place(x = x_offset + 700, y = y_offset -30)
        label_11 = tk.Label(self.win_graph, text = 'Abs retention :', anchor = 'w').place(x = x_offset + 700, y = y_offset)
        label_12 = tk.Label(self.win_graph, text = 'uv reference :', anchor = 'w').place(x = x_offset + 700, y = y_offset + 30)
        label_13 = tk.Label(self.win_graph, text = 'uv absorption :', anchor = 'w').place(x = x_offset + 700, y = y_offset + 60)
        label_14 = tk.Label(self.win_graph, text = 'Gain cross section* :', anchor = 'w').place(x = x_offset + 1000, y = y_offset -30)
        label_15 = tk.Label(self.win_graph, text = 'Gcs * PL_retention :', anchor = 'w').place(x = x_offset + 1000, y = y_offset)

        self.label_0 = tk.Label(self.win_graph, anchor = 'w', textvariable=self.text_filename)
        self.label_0.place(x=x_offset + 120, y = y_offset - 30)
        self.label_1 = tk.Label(self.win_graph, textvariable=self.text_samplename, anchor = 'w')
        self.label_1.place(x = x_offset + 120, y = y_offset)
        self.label_2 = tk.Label(self.win_graph, textvariable=self.text_solvent, anchor = 'w')
        self.label_2.place(x = x_offset + 120, y = y_offset + 30)
        self.label_3 = tk.Label(self.win_graph, textvariable=self.text_absmax, anchor = 'w')
        self.label_3.place(x = x_offset+ 120, y = y_offset + 60)
        self.label_4 = tk.Label(self.win_graph, textvariable=self.text_absend, anchor = 'w')
        self.label_4.place(x = x_offset+ 120, y = y_offset + 90)
        self.label_5 = tk.Label(self.win_graph, textvariable=self.text_PLmax, anchor = 'w')
        self.label_5.place(x = x_offset+ 160 + 350, y = y_offset - 30)
        self.label_6 = tk.Label(self.win_graph, textvariable=self.text_gainmax, anchor = 'w')
        self.label_6.place(x = x_offset+ 160 + 350, y = y_offset)
        self.label_7 = tk.Label(self.win_graph, textvariable=self.text_relativeQY, anchor = 'w')
        self.label_7.place(x = x_offset+ 160 + 350, y = y_offset + 30)
        self.label_8 = tk.Label(self.win_graph, textvariable=self.text_PLlifetime, anchor = 'w')
        self.label_8.place(x = x_offset+ 160 + 350, y = y_offset + 60)
        self.label_9 = tk.Label(self.win_graph, textvariable=self.text_radiative_decay, anchor = 'w')
        self.label_9.place(x = x_offset+ 160 + 350, y = y_offset + 90)
        self.label_10 = tk.Label(self.win_graph, textvariable=self.text_PLretention, anchor = 'w')
        self.label_10.place(x = x_offset+ 120 + 700, y = y_offset -30)
        self.label_11 = tk.Label(self.win_graph, textvariable=self.text_Absretention, anchor = 'w')
        self.label_11.place(x = x_offset+ 120 + 700, y = y_offset)
        self.label_12 = tk.Label(self.win_graph, textvariable=self.text_uvref, anchor = 'w')
        self.label_12.place(x = x_offset+ 120 + 700, y = y_offset + 30)
        self.label_13 = tk.Label(self.win_graph, textvariable=self.text_uvabs, anchor = 'w')
        self.label_13.place(x = x_offset+ 120 + 700, y = y_offset + 60)
        self.label_14 = tk.Label(self.win_graph, textvariable=self.text_gain_cs, anchor = 'w')
        self.label_14.place(x = x_offset+ 140 + 1000, y = y_offset -30)
        self.label_15 = tk.Label(self.win_graph, textvariable=self.text_gcs_PLretention, anchor = 'w')
        self.label_15.place(x = x_offset+ 140 + 1000, y = y_offset)

        self.label_20 = tk.Label(self.win_graph, text = '*Gain cross section was calculated with nF = %s' %nF , anchor = 'w')
        self.label_20.place(x = x_offset+ 1150, y = y_offset + 110)


    def show_job_info(self):

        i = self.r2_val.get()
        d = self.load_data(self.data_list[i]['filename'])

        if 'job' not in d.keys():
            messagebox.showwarning('Warning', 'Job info is not recorded in the file')
        else :
            self.win_job = tk.Toplevel()
            self.win_job.wm_geometry('1000x300')
            self.win_job.title('Job info')

            for i, (key, val) in enumerate(d['job'].items()):
                label = tk.Label(self.win_job, text = '%s : %s' %(key, val), anchor = 'w').place(x = 10, y = 10 + 25 *i) 

    def show_measure_info(self):

        i = self.r2_val.get()
        d = self.load_data(self.data_list[i]['filename'])

        self.win_mesurement = tk.Toplevel()
        self.win_mesurement.wm_geometry('1200x600')
        self.win_mesurement.title('Measurement info')

        i = 0
        for key, val in d['metadata'].items():
            if key not in ['data_path', 'data_processing', 'solvent', 'experiment_time']:
                label = tk.Label(self.win_mesurement, text = '%s :' %key, anchor = 'w').place(x = 20 + 220 * i, y = 50) 
                for j, (key1, val1) in enumerate(val.items()):
                    label = tk.Label(self.win_mesurement, text = '%s : %s' %(key1, val1), anchor = 'w').place(x = 40 + 220 * i, y = 75 + 25 *j) 
                i += 1
            elif key == 'experiment_time':
                label = tk.Label(self.win_mesurement, text = '%s : %s' %(key, val), anchor = 'w').place(x = 20, y = 10) 
            elif key == 'data_processing':
                label = tk.Label(self.win_mesurement, text = '%s :' %key, anchor = 'w').place(x = 20, y = 425) 
                for j, (key1, val1) in enumerate(val.items()):
                    label = tk.Label(self.win_mesurement, text = '%s : %s' %(key1, val1), anchor = 'w').place(x = 40, y = 450 + 25 *j)     

    def show_heatmap(self):

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        from matplotlib.colors import Normalize 

        self.win_heatmap = tk.Toplevel()
        self.win_heatmap.wm_geometry('1000x600')
        self.win_heatmap.title('Heat Map')

        data_sorted = sorted(self.data_list, key=lambda x:x['PL_lambda_max'])        
        data_list = []

        for data in data_sorted:
            d = []
            for key, val in self.config.items():
                if data[key]:
                    d.append(data[key]/val['Normalize_at'])
                else:
                    d.append(0)
            data_list.append(d)

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(1,1,1)

        # color = "gist_rainbow_r"
        color = "summer"
        # color = "nipy_spectral"
        # color = "gist_heat"

        im, cbar = heatmap(np.asarray(data_list), [d['name'] for d in data_sorted], self.spoke_labels, ax=ax,
                           cbar_kw={'shrink' :0.75}, cmap=color, norm=Normalize(vmin=0, vmax=1), aspect = 4/len(data_list),  cbarlabel=None)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, master= self.win_heatmap)
        canvas.draw()
        canvas.get_tk_widget().place(x=10, y = 10, height=650, width=950)


    def show_scatter_plots(self, label):

        fig = plt.figure(figsize=(13,7))

        ax = []
        ylabels = ['spectal factor / cm2 s', 'PL quantum yield', 'radiative decay constant / s-1', 'Gain cross section / cm2', 'PL retention', 'Gain cross section * PL retention / cm2']
        titles = ['Spectral factor', 'PL quantum yield', 'Radiative decay constant', 'Gain cross section', 'PL retention', 'Gain cross section\n* PL retention']
        data = ['max_gain_factor', 'relative_QY', 'radiative_decay', 'Gain_cross_section', 'PL_maintenance', 'Gain*PL_retention' ]

        for i in range(6):
            ax.append(fig.add_subplot(2,3,i + 1, xlabel = 'lambda_max(PL)/ nm', ylabel = ylabels[i], xlim = (400,700)))
            ax[i].set_title(titles[i])
            for d in self.data_list:
                ax[i].scatter(d['PL_lambda_max'], d[data[i]])
                if label:
                    ax[i].text(d['PL_lambda_max'], d[data[i]], d['name'], fontsize = 5)

        fig.tight_layout()
        fig.show()



    def calc_color(self, PL_spectrum):        
        
        SC = Spec_to_Color()
        rgb = SC.spec_to_rgb(PL_spectrum)
        color = SC.hex_to_str(rgb)
        return color


    def update_chart(self, data_num):

        data = self.data_list[data_num]

        self.fig2.clear()
        self.ax7 = self.fig2.add_subplot(1,1,1, projection = 'radar')

        d = []

        for key, val in self.config.items():
            if data[key]:
                d.append(data[key]/val['Normalize_at'])
            else:
                d.append(0)

        if data['PL_rgb']:
            color = data['PL_rgb']
        else : color = 'b'

        self.ax7.plot(self.theta, d, color=color)
        self.ax7.fill(self.theta, d, facecolor=color, alpha=0.25)

        self.ax7.set_rgrids([0.2, 0.4, 0.6, 0.8, 1.0], labels = ['','','','',''])
        self.ax7.set_varlabels(self.spoke_labels)

        self.canvas2.draw()


    def update_plot(self, data_num):

        d = self.load_data(self.data_list[data_num]['filename'])

        cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

        self.ax1.clear()
        self.ax2.clear()
        self.ax3.clear()
        self.ax4.clear()
        self.ax5.clear()
        self.ax6.clear()
        self.ax1.set_xlabel('wavelength / nm')
        self.ax1.set_ylabel('Absorbance')
        self.ax2.set_ylabel('PL emission/energy nm-1 s-1')
        self.ax3.set_xlabel('wavelength / nm')
        self.ax3.set_ylabel('Absorbance')
        self.ax4.set_ylabel('Gain factor/1E24 cm2 s')
        self.ax5.set_xlabel('time / ns')
        self.ax5.set_ylabel('photon count')
        self.ax5.set_yscale('log')
        self.ax6.set_xlabel('time / s')
        self.ax6.set_ylabel('Normalized Absorbance(365nm)/PL intensity')

        if d['absorption'] is not None:

            calc_range = self._to_index_range(d['absorption']['absorbance'][0], d['metadata']['absorption']['abs_calc_range'])

            self.ax1.plot(d['absorption']['absorbance'][0], d['absorption']['absorbance'][1], color = cycle[0])
            self.ax1.set_ylim((-0.2* np.max(d['absorption']['absorbance'][1][calc_range[0]:calc_range[1]+1]), 1.2 * np.max(d['absorption']['absorbance'][1][calc_range[0]:calc_range[1]+1])))
            self.ax3.plot(d['absorption']['absorbance'][0], d['absorption']['absorbance'][1], color = cycle[0])
            self.ax3.set_ylim((-0.2* np.max(d['absorption']['absorbance'][1][calc_range[0]:calc_range[1]+1]), 1.2 * np.max(d['absorption']['absorbance'][1][calc_range[0]:calc_range[1]+1])))
        
        if d['PL'] is not None:

            calc_range =  self._to_index_range(d['PL']['energy'][0], d['metadata']['PL']['PL_calc_range'])

            self.ax2.plot(d['PL']['energy'][0], d['PL']['energy'][1], color = cycle[1])
            self.ax2.set_xlim((300,800))
            self.ax2.set_ylim((-0.2* np.max(d['PL']['energy'][1][calc_range[0]:calc_range[1]]), 1.2 * np.max(d['PL']['energy'][1][calc_range[0]:calc_range[1]])))
            self.ax4.plot(d['PL']['gain_spectrum'][0], d['PL']['gain_spectrum'][1]*1E24, color = cycle[1])
            self.ax4.set_xlim((300,800))
            self.ax4.set_ylim((-0.2* np.max(d['PL']['gain_spectrum'][1][calc_range[0]:calc_range[1]]*1E24), 1.2 * np.max(d['PL']['gain_spectrum'][1][calc_range[0]:calc_range[1]]*1E24)))
            if d['uv']['absorbance_time_trace'] is not None:
                self.ax6.plot(d['uv']['absorbance_time_trace'][0], d['uv']['absorbance_time_trace'][1]/d['uv']['absorbance_time_trace'][1,0], label = 'Abs', color = cycle[0])
            try:
                try:
                    self.ax6.plot(d['PL']['int_time_trace'][0], d['PL']['int_time_trace'][1]/d['PL']['int_time_trace'][1,0], label = 'PL', color = cycle[1])
                except:
                    self.ax6.plot(d['PL']['peak_time_trace'][0], d['PL']['peak_time_trace'][1]/d['PL']['peak_time_trace'][1,0], label = 'PL', color = cycle[1])
            except: pass
            self.ax6.legend()

        if d['TE'] is not None:
            self.ax5.plot(d['TE']['raw_data'][0], d['TE']['raw_data'][1], color = cycle[0])
            self.ax5.plot(d['TE']['fitted_data'][0], d['TE']['fitted_data'][1], color = cycle[1])
            try:
                self.ax5.set_xlim((10, 0.9*1E9/d['TE']['metadata']['excitation_frequency']))
            except:
                self.ax5.set_xlim((10, 0.9*1E9/d['metadata']['TE']['excitation_frequency']))
            self.ax5.set_ylim((1E0,np.max(d['TE']['raw_data'][1])*2))

        self.canvas1.draw()



        


    def show_mol_image(self, smiles):

        global img

        mol = Chem.MolFromSmiles(smiles)
        tm = rdMolDraw2D.PrepareMolForDrawing(mol)

        d = rdMolDraw2D.MolDraw2DCairo(300, 250)
        d.drawOptions().minFontSize = 13
        d.SetFontSize(1.2)
        d.DrawMolecule(tm)     
        d.FinishDrawing()
        d.WriteDrawingText(image_file) 

        img = tk.PhotoImage(file=image_file)
        self.canvas0.create_image(10, 15, image=img, anchor=tk.NW)


    def select_data(self):

        i = self.r2_val.get()
        d = self.data_list[i]

        #update window
        self.text_filename.set(os.path.basename(d['filename']))
        self.text_samplename.set(d['name'])
        self.text_solvent.set(d['solvent'])

        self.text_absmax.set('{:.4f} ({:.1f} nm)'.format(d['abs_max'], d['abs_lambda_max']))
        self.text_absend.set('{:.1f} nm'.format(d['abs_lambda_end']))

        self.text_uvref.set('{:.3f} mW'.format(d['uv_reference']*1000))
        self.text_uvabs.set('{:.3f} mW'.format(d['uv_absorption']*1000))
        self.text_PLmax.set('{:.2f} ({:.1f} nm)'.format(d['PL_max'], d['PL_lambda_max']))
        self.text_relativeQY.set('{:.3f}'.format(d['relative_QY'])) ###
        self.text_gainmax.set('{:.3f} ({:.2f} nm)'.format(d['max_gain_factor']*1E24, d['max_gain_wavelength']))

        if d['PL_maintenance']:
            self.text_PLretention.set('{:.3f}'.format(d['PL_maintenance']))
        else : 
            self.text_PLretention.set(d['PL_maintenance'])

        if d['uv_maintenance']:
            self.text_Absretention.set('{:.3f}'.format(d['uv_maintenance']))
        else : 
            self.text_Absretention.set(d['uv_maintenance'])
        
        if d['PL_lifetime']:
            self.text_PLlifetime.set('{:.3f} ns'.format(d['PL_lifetime']))
        else:
            self.text_PLlifetime.set(d['PL_lifetime'])

        if d['radiative_decay']:
            self.text_radiative_decay.set('{:.3e} s-1'.format(d['radiative_decay']))
        else:
            self.text_radiative_decay.set(d['radiative_decay'])

        if d['Gain_cross_section']:
            self.text_gain_cs.set('{:.3e} cm2'.format(d['Gain_cross_section']))
        else:
            self.text_gain_cs.set(d['Gain_cross_section'])

        if d['Gain*PL_retention']:
            self.text_gcs_PLretention.set('{:.3e} cm2'.format(d['Gain*PL_retention']))
        else:
            self.text_gcs_PLretention.set(d['Gain*PL_retention'])

        if d['smiles']:
            self.show_mol_image(d['smiles'])
    
        self.update_plot(i)
        self.update_chart(i)


    def load_data(self, filename):
        with open(filename, 'rb') as f:
            d = pkl.load(f)
        return d


    def _to_index_range(self, wavelength, wavelength_range):

        return [np.where(wavelength >= wavelength_range[0])[0][0], np.where(wavelength >= wavelength_range[1])[0][0]]


    def get_data(self):

        self.data_list = []
        #filename
        for filename in self.filenames:
            d = self.load_data(filename)

            keys = ['filename', 'name', 'solvent', 'smiles', 'abs_max', 'abs_lambda_max', 'abs_lambda_end', 'uv_reference', 'uv_absorption', 'uv_maintenance', \
                    'PL_max', 'PL_lambda_max', 'PL_maintenance', 'PL_rgb', 'relative_QY', 'max_gain_factor', 'max_gain_wavelength', 'PL_lifetime', \
                    'radiative_decay', 'Gain_cross_section', 'Gain*PL_retention', 'abs_calc_range', 'PL_calc_range']
            data = {key : None for key in keys}

            data['filename'] = filename

        #sample name
            if 'job' in d.keys():
                data['name'] = d['job']['target_name']
            elif 'sample' in d['metadata'].keys():
                data['name'] = d['metadata']['sample']['name']
        #solvent
            if 'solvent' in d['metadata'].keys():
                data['solvent'] = d['metadata']['solvent']
            elif 'sample' in d['metadata'].keys():
                data['solvent'] = d['metadata']['sample']['solvent']
        #smiles
            if 'job' in d.keys() and 'target_smiles' in d['job'].keys():
                data['smiles'] = d['job']['target_smiles']
        #data
            if d['absorption'] is not None:
                data['abs_max'] = d['absorption']['max']
                data['abs_lambda_max'] = d['absorption']['lambda_max']
                data['abs_lambda_end'] = d['absorption']['end_wavelength']
            if d['PL'] is not None:
                data['uv_reference'] = d['uv']['reference']
                data['uv_absorption'] = d['uv']['absorption']
                data['PL_max'] = d['PL']['max']
                data['PL_lambda_max'] = d['PL']['lambda_max']
                data['relative_QY'] = d['PL']['relative_QY']
                if d['uv']['absorbance_maintenance(at_1min)'] is not None:
                    data['uv_maintenance'] = d['uv']['absorbance_maintenance(at_1min)']
                if d['PL']['maintenance(at_1min)'] is not None:
                    data['PL_maintenance'] = d['PL']['maintenance(at_1min)']
                if  d['PL']['max_gain_factor'] is not None:
                    data['max_gain_factor'] = d['PL']['max_gain_factor']
                    data['max_gain_wavelength'] = d['PL']['max_gain_wavelength']
                data['PL_rgb'] = self.calc_color(d['PL']['energy'])
            if d['TE'] is not None:
                if 'tau' in d['TE']['fitting_results'].keys():
                    data['PL_lifetime'] = d['TE']['fitting_results']['tau']
                elif 'tau1' in d['TE']['fitting_results'].keys():
                    if d['TE']['fitting_results']['amp1'] > d['TE']['fitting_results']['amp2']:
                        data['PL_lifetime'] = d['TE']['fitting_results']['tau1']
                    else:
                        data['PL_lifetime'] = d['TE']['fitting_results']['tau2']
            if data['PL_lifetime'] and data['relative_QY']:
                data['radiative_decay'] = 1E9 * data['relative_QY']/data['PL_lifetime']
            if data['radiative_decay'] and data['max_gain_factor']:
                data['Gain_cross_section'] = data['radiative_decay'] * data['max_gain_factor'] /nF**2
            if data['Gain_cross_section'] and data['PL_maintenance']:
                data['Gain*PL_retention'] = data['Gain_cross_section'] * data['PL_maintenance']


            self.data_list.append(data) 


    def file_open(self):
        self.filenames = filedialog.askopenfilenames(filetypes=[("pickle file", ".pkl")])

        self.get_data()

        self.open_graph_window()

        if self.r_buttons:
            for i in range(len(self.r_buttons)):
                self.r_buttons[-i-1].destroy()

        self.r_buttons = []
        for i, filename in enumerate(self.filenames):
            fname = os.path.basename(filename)
            row = i%36
            column = i//36
            self.r_buttons.append(tk.Radiobutton(self.win_graph, text = fname.split('.')[0], value = i, var = self.r2_val, command = self.select_data))
            self.r_buttons[i].place(x = 1350 + 150 * column, y = 40 + 22 * row)

        self.select_data()

    def exit(self):
        self.master.quit()
        
root = tk.Tk()
main = Application(master=root)
main.mainloop()

 

    # def get_sample_names(self):

    #     sample_names = []

    #     for filename in self.filenames:
    #         with open(filename, 'rb') as f:
    #             d = pkl.load(f)

    #         if 'job' in d.keys():
    #             sample_names.append(d['job']['target_name'])
    #         elif 'sample' in d['metadata'].keys():
    #             sample_names.append(d['metadata']['sample']['name'])
    #         else :
    #             sample_names.append('No name')
    #     return sample_names

    # def load_data(self, filename):

    #     with open(filename, 'rb') as f:
    #         self.data = pkl.load(f)
    #     #shoe filename
        


    #     # self.show_mol_image('C1(N2C(C=CC=C3)=C3C4=C2C=CC=C4)=CC=C(/C=C/C5=CC=C(C6=CC=C(/C=C/C7=CC=C(N8C(C=CC=C9)=C9C%10=C8C=CC=C%10)C=C7)C=C6)C=C5)C=C1')
    #     # self.canvas0.configure(background = self.calc_color())
    #     self.update_plot()
    #     self.update_chart()



    # def show_scatter_plots(self):

    #     from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    #     def on_closing():
    #         plt.close()
    #         self.win_scatter.destroy()

    #     def update_scatter():

    #         text = c_val.get()

    #         for i in range(6):
    #             ax[i].clear()
    #             ax[i].set_xlabel('lambda_max(PL)/ nm')
    #             ax[i].set_ylabel(ylabels[i])
    #             ax[i].set_title(titles[i])
    #             ax[i].set_xlim((400,700))

    #         for d in self.data_list:
    #             for i in range(6):
    #                 ax[i].scatter(d['PL_lambda_max'], d[data[i]])
    #                 if text:
    #                     ax[i].text(d['PL_lambda_max'], d[data[i]], d['name'], fontsize = 5)

    #         canvas.draw()

    #     self.win_scatter = tk.Toplevel()
    #     self.win_scatter.wm_geometry('1200x800')
    #     self.win_scatter.title('Scatter plots')

    #     c_val = tk.BooleanVar()
        
    #     c = tk.Checkbutton(self.win_scatter, variable = c_val, text = 'sample_number', command = update_scatter)
    #     c.place(x = 1000, y = 10)

    #     fig = plt.figure(figsize=(13,7))

    #     ax = []
    #     ylabels = ['spectal factor / cm2 s', 'PL quantum yield', 'radiative decay constant / s-1', 'Gain cross section / cm2', 'PL retention', 'Gain cross section * PL retention / cm2']
    #     titles = ['Spectral factor', 'PL quantum yield', 'Radiative decay constant', 'Gain cross section', 'PL retention', 'Gain cross section\n* PL retention']
    #     data = ['max_gain_factor', 'relative_QY', 'radiative_decay', 'Gain_cross_section', 'PL_maintenance', 'Gain*PL_retention' ]

    #     for i in range(6):
    #         ax.append(fig.add_subplot(2,3,i + 1, xlabel = 'lambda_max(PL)/ nm', ylabel = ylabels[i], xlim = (400,700)))
    #         ax[i].set_title(titles[i])

    #     fig.tight_layout()

    #     canvas = FigureCanvasTkAgg(fig, master= self.win_scatter)
    #     canvas.get_tk_widget().place(x=10, y = 30, height=700, width=1100) 

    #     update_scatter()

    #     self.win_scatter("WM_DELETE_WINDOW", on_closing) 