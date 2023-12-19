import numpy as np
import matplotlib.pyplot as plt
import json
import pickle as pkl
import os
from tkinter import filedialog
import tkinter as tk
from tkinter import messagebox, ttk
import copy

class Application(tk.Frame):

    def __init__(self, master = None):
        super().__init__(master, height = 850, width = 1500)
        self.pack()

        #self.master.geometry("300x300")
        self.master.title("Optical Characterization results")

        ### Menu
        menu_bar = tk.Menu(self, tearoff=0)
        # File
        menu_file = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File", menu=menu_file,  underline=0)
        menu_file.add_command(label="Open", command=self.file_open, underline=0, accelerator = 'Ctrl-O')
        menu_file.add_command(label="Exit", command=self.exit, underline=0, accelerator = 'Ctrl-Q')
        # File info
        menu_info = tk.Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label="File_info", menu=menu_info,  underline=0)
        # menu_info.add_command(label="Injection_info", command=self.show_injection_info, underline=0)
        # menu_info.add_command(label="HPLC_info", command=self.show_HPLC_info, underline=0)
        menu_info.add_command(label="Analysis_info", command=self.show_analysis_info, underline=0)

        # menu_plots = tk.Menu(menu_bar, tearoff=0)
        # menu_bar.add_cascade(label="Overlay", menu=menu_plots,  underline=0)
        # menu_plots.add_command(label="chromatogram", command=self.DAD_plots, underline=0)


        self.master.config(menu=menu_bar)
        self.r2_val = tk.IntVar()
        self.text_filename = tk.StringVar()
        self.text_samplename = tk.StringVar()
        self.text_absmax = tk.StringVar()
        self.text_absend = tk.StringVar()
        self.text_PLmax = tk.StringVar()
        self.text_gainmax = tk.StringVar()
        self.text_relativeQY = tk.StringVar()
        self.text_PLlifetime = tk.StringVar()
        self.text_PLretention = tk.StringVar()
        self.text_Absretention = tk.StringVar()
        self.text_uvref  = tk.StringVar()
        self.text_uvabs = tk.StringVar()

        self.data = None
        self.r_buttons = None
        

        self.create_widgets()

        self.master.protocol("WM_DELETE_WINDOW", self.exit)

    #wrapper
    # def data_method(func):
    #     def wrapper(self, *args, **kwargs):
    #         if not self.data:
    #             messagebox.showwarning('Warning', 'Please load results file')
    #         else: 
    #             func(self, *args, **kwargs)
    #     return wrapper   
    

    def create_widgets(self):

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        #filename
        self.label_0 = tk.Label(self, text = 'filename :', anchor = 'w')
        self.label_0.place(x=10, y = 10, height=20, width=100)

        ###PL transient
        # self.fig3 = plt.figure()
        # self.ax5 = self.fig3.add_subplot(1,1,1, xlabel = 'time / ns', ylabel = 'photon count')
        # self.canvas3 = FigureCanvasTkAgg(self.fig3, master= self)
        # self.canvas3.draw()
        # self.canvas3.get_tk_widget().place(x=30, y = 400, height=380, width=500)

        # ###Abs/PL time trace
        # self.fig4 = plt.figure()
        # self.ax6 = self.fig4.add_subplot(1,1,1, xlabel = 'time / s', ylabel = 'Normalized Absorbance/PL intensity')
        # self.canvas4 = FigureCanvasTkAgg(self.fig4, master= self)
        # self.canvas4.draw()
        # self.canvas4.get_tk_widget().place(x=540, y = 400, height=380, width=500)

        ###Absorption and PL spectra
        self.fig1 = plt.figure(figsize=(8,8))
        self.ax1 = self.fig1.add_subplot(2,2,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance') #Absorption/PL 
        self.ax2 = self.ax1.twinx()
        self.ax2.set_ylabel('PL emission/energy nm-1 s-1')
        self.ax3 = self.fig1.add_subplot(2,2,2, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
        self.ax4 = self.ax3.twinx()
        self.ax4.set_ylabel('Gain factor/1E24 cm2 s')
        self.ax5 = self.fig1.add_subplot(2,2,3, xlabel = 'time / ns', ylabel = 'photon count')
        self.ax6 = self.fig1.add_subplot(2,2,4, xlabel = 'time / s', ylabel = 'Normalized Absorbance(365nm)/PL intensity')

        self.fig1.tight_layout()

        self.canvas1 = FigureCanvasTkAgg(self.fig1, master= self)
        self.canvas1.draw()
        # self.canvas1.get_tk_widget().place(x=30, y = 30, height=380, width=420)
        self.canvas1.get_tk_widget().place(x=30, y = 30, height=600, width=1000)

        x_offset = 1050
        y_offset = 50
        label_1 = tk.Label(self, text = 'sample name :', anchor = 'w').place(x = x_offset, y = y_offset)
        label_2 = tk.Label(self, text = 'Absorption max :', anchor = 'w').place(x = x_offset, y = y_offset + 30)
        label_3 = tk.Label(self, text = 'Absorption end :', anchor = 'w').place(x = x_offset, y = y_offset + 60)
        label_4 = tk.Label(self, text = 'PL max :', anchor = 'w').place(x = x_offset, y = y_offset + 90)
        label_5 = tk.Label(self, text = 'Gain factor max(1E-24cm2 s) :', anchor = 'w').place(x = x_offset, y = y_offset + 120)
        label_6 = tk.Label(self, text = 'Relative QY :', anchor = 'w').place(x = x_offset, y = y_offset + 150)
        label_7 = tk.Label(self, text = 'PL lifetime (ns) :', anchor = 'w').place(x = x_offset, y = y_offset + 180)
        label_8 = tk.Label(self, text = 'PL retention :', anchor = 'w').place(x = x_offset, y = y_offset + 210)
        label_9 = tk.Label(self, text = 'Abs retention :', anchor = 'w').place(x = x_offset, y = y_offset + 240)
        label_10 = tk.Label(self, text = 'uv reference (mW) :', anchor = 'w').place(x = x_offset, y = y_offset + 270)
        label_11 = tk.Label(self, text = 'uv absorption (mW) :', anchor = 'w').place(x = x_offset, y = y_offset + 300)

        self.label_1 = tk.Label(self, textvariable=self.text_samplename, anchor = 'w')
        self.label_1.place(x = x_offset + 150, y = y_offset)
        self.label_2 = tk.Label(self, textvariable=self.text_absmax, anchor = 'w')
        self.label_2.place(x = x_offset+ 150, y = y_offset + 30)
        self.label_3 = tk.Label(self, textvariable=self.text_absend, anchor = 'w')
        self.label_3.place(x = x_offset+ 150, y = y_offset + 60)
        self.label_4 = tk.Label(self, textvariable=self.text_PLmax, anchor = 'w')
        self.label_4.place(x = x_offset+ 150, y = y_offset + 90)
        self.label_5 = tk.Label(self, textvariable=self.text_gainmax, anchor = 'w')
        self.label_5.place(x = x_offset+ 150, y = y_offset + 120)
        self.label_6 = tk.Label(self, textvariable=self.text_relativeQY, anchor = 'w')
        self.label_6.place(x = x_offset+ 150, y = y_offset + 150)
        self.label_7 = tk.Label(self, textvariable=self.text_PLlifetime, anchor = 'w')
        self.label_7.place(x = x_offset+ 150, y = y_offset + 180)
        self.label_8 = tk.Label(self, textvariable=self.text_PLretention, anchor = 'w')
        self.label_8.place(x = x_offset+ 150, y = y_offset + 210)
        self.label_9 = tk.Label(self, textvariable=self.text_Absretention, anchor = 'w')
        self.label_9.place(x = x_offset+ 150, y = y_offset + 240)
        self.label_10 = tk.Label(self, textvariable=self.text_uvref, anchor = 'w')
        self.label_10.place(x = x_offset+ 150, y = y_offset + 270)
        self.label_11 = tk.Label(self, textvariable=self.text_uvabs, anchor = 'w')
        self.label_11.place(x = x_offset+ 150, y = y_offset + 300)
        ###Absorption and Gain spectra
        # self.fig2 = plt.figure()
        # self.ax3 = self.fig2.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
        # self.ax4 = self.ax3.twinx()
        # self.ax4.set_ylabel('Gain factor/1E24 cm2 s')
        # self.canvas2 = FigureCanvasTkAgg(self.fig2, master= self)
        # self.canvas2.draw()
        # self.canvas2.get_tk_widget().place(x=540, y = 30, height=380, width=500)



        #plot results###########
        # results_x_offset = 840
        # results_y_offset = 400
        # self.label_4 = tk.Label(self, text = 'plot results :', anchor = 'w')
        # self.label_4.place(x=results_x_offset, y = 20 + results_y_offset, height=20, width=100) 

        # self.checkbox_1 = tk.Checkbutton(self, variable = self.c1_val, text = 'raw_chromatogram', command = self.update_plot)
        # self.checkbox_1.place(x = results_x_offset + 10, y = 50 + results_y_offset)

        # self.checkbox_2 = tk.Checkbutton(self, variable = self.c2_val, text = 'peak number', command = self.update_plot)
        # self.checkbox_2.place(x = results_x_offset + 10, y = 70 + results_y_offset)

        # self.checkbox_3 = tk.Checkbutton(self, variable = self.c3_val, text = 'peak area', command = self.update_plot)
        # self.checkbox_3.place(x = results_x_offset + 10, y = 90 + results_y_offset)



        #plot spectra##########
        # spectra_x_offset = 840
        # spectra_y_offset = 550
        # self.label_5 = tk.Label(self, text = 'plot spectra :', anchor = 'w')
        # self.label_5.place(x=spectra_x_offset, y = 10 + spectra_y_offset, height=20, width=100)

        # r1_text = ['Peak No.', 'Time', 'All']
        # r1 = []
        # for i in range(3):
        #     r1.append(tk.Radiobutton(self, text=r1_text[i], value=i, var=self.r1_val))
        #     r1[i].place(x = spectra_x_offset + 20, y = 60 + 25 * (i-1) + spectra_y_offset)

        # self.entry_peak_num = tk.Entry(self, textvariable = self.peak_num, width = 6)
        # self.entry_peak_num.place(x = spectra_x_offset + 120, y = 35 + spectra_y_offset)
        # self.entry_r_time = tk.Entry(self, textvariable = self.retention_time, width = 6)
        # self.entry_r_time.place(x = spectra_x_offset + 120, y = 60 + spectra_y_offset)

        # self.button_DAD_spec = tk.Button(self,  text = 'plot', command = self.plot_absorption, height = 2, width = 7)
        # self.button_DAD_spec.place(x = spectra_x_offset + 170, y = 35 + spectra_y_offset)


    def DAD_plots(self):

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        def on_closing():
            plt.close()
            self.win_DAD_chromatograms.destroy()

        def update_plot():
            if r1_val.get() == 0:
                chromatograms = r_chromatograms
            else:
                chromatograms = c_chromatograms 

            ax1.clear()
            ax1.set_xlabel('time / min')
            ax1.set_ylabel('mAU')
            for i in range(len(self.filenames)):
                if c_vals[i].get() == True:
                    ax1.plot(chromatograms[i][0], chromatograms[i][1], label = injection_names[i])
            ax1.legend()
            canvas.draw()
        
        def clear_select():
            for i in range(len(self.filenames)):
                c_vals[i].set(False)
            ax1.clear()
            canvas.draw()

        r_chromatograms = []
        c_chromatograms = []
        injection_names = []

        fig = plt.figure(figsize = (14, 4))
        ax1 = fig.add_subplot(1,1,1, xlabel = 'time / min', ylabel = 'mAU')

        for filename in self.filenames:

            with open(filename, 'rb') as f:
                d = pkl.load(f)

            r_chromatograms.append(d['chromatogram'])
            c_chromatograms.append(d['corrected_chromatogram'])
            injection_names.append(os.path.basename(filename).split('.')[0])

        self.win_DAD_chromatograms = tk.Toplevel()
        self.win_DAD_chromatograms.wm_geometry('1700x500')
        self.win_DAD_chromatograms.title('DAD_chromatograms')

        r1_text = ['raw', 'bkg_corrected']
        r1 = []
        r1_val = tk.IntVar()
        for i in range(2):
            r1.append(tk.Radiobutton(self.win_DAD_chromatograms, text=r1_text[i], value=i, var=r1_val, command = update_plot))
            r1[i].place(x = 1150 + i * 50, y = 10)

        checkboxs = []
        c_vals = []
        for i in range(len(self.filenames)):
            row = i%15
            column = i//15
            c_vals.append(tk.BooleanVar())
            checkboxs.append(tk.Checkbutton(self.win_DAD_chromatograms, variable = c_vals[i], text = injection_names[i], command = update_plot))
            checkboxs[i].place(x = 1400 + 120 * column, y = 50 + 22 * row)

        button_clear_select = tk.Button(self.win_DAD_chromatograms,  text = 'clear', command = clear_select, height = 2, width = 7)
        button_clear_select.place(x = 1400, y = 385)

        canvas = FigureCanvasTkAgg(fig, master= self.win_DAD_chromatograms)
        canvas.draw()
        canvas.get_tk_widget().place(x = 0, y = 30)

        self.win_DAD_chromatograms.protocol("WM_DELETE_WINDOW", on_closing)


    @data_method
    def show_analysis_info(self):

        if 'analysis_setting' not in self.data.keys():
            messagebox.showwarning('warning', 'analysis info is not stored in the file')
            
        else:
            self.win_Anal_info = tk.Toplevel()
            self.win_Anal_info.wm_geometry('400x350')
            self.win_Anal_info.title('analysis setting')
            for i, (key, val) in enumerate(self.data['analysis_setting'].items()):
                if 'plot' not in key and 'varbose' not in key:
                    text = '%s : %s' %(key, val) 
                    label = tk.Label(self.win_Anal_info, text = text, anchor = 'w')
                    label.place(x = 10, y = 10 + 25*i)


    @data_method
    def plot_DAD_chromatogram(self):
        wavelength = self.wavelength.get()
        chromatogram = self.Chrom.get_DAD_chromatogram(self.data['3D_fields'], wavelength - 1, wavelength + 1)

        fig = plt.figure()
        ax = fig.add_subplot(1,1,1, title = 'DAD chromatogram at %s nm' %wavelength, xlabel = 'time/min', ylabel = 'absorbance/mAU')
        ax.plot(chromatogram[0], chromatogram[1])
        fig.show()


    @data_method
    def plot_absorption(self):
        mode = self.r1_val.get()
        if mode == 0:
            peak = self.peak_num.get()
            peak_data = self.data['peak_data'][peak]

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1, title = 'Absorption spectrum at peak %s' %peak, xlabel = 'wavelength/nm', ylabel = 'absorbance/mAU')
            ax.plot(peak_data['spectrum'][0], peak_data['spectrum'][1])
            ax.legend(['{:.3f} min'.format(peak_data['peak']['time'])], loc = 'upper right')
            fig.show()
        if mode == 1:
            r_time = self.retention_time.get()
            spectrum = self.Chrom.get_DAD_spectrum(self.data['3D_fields'], retention_time=r_time)

            fig = plt.figure()
            ax = fig.add_subplot(1,1,1, title = 'Absorption spectrum at %s min' %r_time, xlabel = 'wavelength/nm', ylabel = 'absorbance/mAU')
            ax.plot(spectrum[0], spectrum[1])
            fig.show()
        if mode == 2:
            self.plot_peak_spectra(self.data['peak_data'])

    @data_method
    def plot_peak_spectra(self, peak_data, y_label = 'mAU'):

        import math
        from matplotlib import gridspec 

        row_num = math.ceil(len(peak_data) / 5)
        figs1 =plt.figure(figsize = (20, row_num * 2.5))
        gs = gridspec.GridSpec(row_num ,5)

        ax = [[] for i in range(len(peak_data))]
        for i, peak in enumerate(peak_data):
            ax[i] = figs1.add_subplot(gs[i //5, i % 5], xlabel = 'wavelength/nm', ylabel = y_label)
            ax[i].plot(peak['spectrum'][0], peak['spectrum'][1])
            ax[i].legend(['{:.3f} min'.format(peak['peak']['time'])], loc = 'upper right')
        figs1.tight_layout() 
        figs1.show()


    # self.fig3 = plt.figure()
    #     self.ax5 = self.fig3.add_subplot(1,1,1, xlabel = 'time / ns', ylabel = 'photon count')
    #     self.canvas3 = FigureCanvasTkAgg(self.fig3, master= self)
    #     self.canvas3.draw()
    #     self.canvas3.get_tk_widget().place(x=30, y = 400, height=380, width=480)

    #     ###Abs/PL time trace
    #     self.fig4 = plt.figure()
    #     self.ax6 = self.fig4.add_subplot(1,1,1, xlabel = 'time / s', ylabel = 'Normalized Absorbance/PL intensity')
    #     self.canvas4 = FigureCanvasTkAgg(self.fig4, master= self)
    #     self.canvas4.draw()
    #     self.canvas4.get_tk_widget().place(x=520, y = 400, height=380, width=480)

    #     ###Absorption and PL spectra
    #     self.ax2.set_ylabel('PL emission/energy nm-1 s-1')
    #     # self.ax3 = self.fig1.add_subplot(2,2,2, xlabel = 'wavelength', ylabel = 'mAU') #Gain/PL
    #     # self.fig1.tight_layout()
    #     self.canvas1 = FigureCanvasTkAgg(self.fig1, master= self)
    #     self.canvas1.draw()
    #     # self.canvas1.get_tk_widget().place(x=30, y = 30, height=380, width=420)
    #     self.canvas1.get_tk_widget().place(x=30, y = 30, height=380, width=480)

    #     ###Absorption and Gain spectra
    #     self.fig2 = plt.figure()
    #     self.ax3 = self.fig2.add_subplot(1,1,1, xlabel = 'wavelength / nm', ylabel = 'Absorbance')
    #     self.ax4 = self.ax3.twinx()
    #     self.ax4.set_ylabel('Gain factor/1E24 cm2 s')
    #     self.canvas2 = FigureCanvasTkAgg(self.fig2, master= self)
    #     self.canvas2.draw()
    #     self.canvas2.get_tk_widget().place(x=520, y = 30, height=380, width=480)

    def update_plot(self):

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

        if self.data['absorption'] is not None:
            self.ax1.plot(self.data['absorption']['absorbance'][0], self.data['absorption']['absorbance'][1], color = cycle[0])
            self.ax1.set_ylim((-0.2* np.max(self.data['absorption']['absorbance'][1]), 1.2 * np.max(self.data['absorption']['absorbance'][1])))
            self.ax3.plot(self.data['absorption']['absorbance'][0], self.data['absorption']['absorbance'][1], color = cycle[0])
            self.ax3.set_ylim((-0.2* np.max(self.data['absorption']['absorbance'][1]), 1.2 * np.max(self.data['absorption']['absorbance'][1])))
        
        if self.data['PL'] is not None:
            self.ax2.plot(self.data['PL']['energy'][0], self.data['PL']['energy'][1], color = cycle[1])
            self.ax2.set_xlim((300,800))
            self.ax2.set_ylim((-0.2* np.max(self.data['PL']['energy'][1]), 1.2 * np.max(self.data['PL']['energy'][1])))
            self.ax4.plot(self.data['PL']['gain_spectrum'][0], self.data['PL']['gain_spectrum'][1]*1E24, color = cycle[1])
            self.ax4.set_xlim((300,800))
            self.ax4.set_ylim((-0.2* np.max(self.data['PL']['gain_spectrum'][1]*1E24), 1.2 * np.max(self.data['PL']['gain_spectrum'][1]*1E24)))
            if self.data['uv']['absorbance_time_trace'] is not None:
                self.ax6.plot(self.data['uv']['absorbance_time_trace'][0], self.data['uv']['absorbance_time_trace'][1]/self.data['uv']['absorbance_time_trace'][1,0], label = 'Abs', color = cycle[0])
            try:
                try:
                    self.ax6.plot(self.data['PL']['int_time_trace'][0], self.data['PL']['int_time_trace'][1]/self.data['PL']['int_time_trace'][1,0], label = 'PL', color = cycle[1])
                except:
                    self.ax6.plot(self.data['PL']['peak_time_trace'][0], self.data['PL']['peak_time_trace'][1]/self.data['PL']['peak_time_trace'][1,0], label = 'PL', color = cycle[1])
            except: pass
            self.ax6.legend()

        if self.data['TE'] is not None:
            self.ax5.plot(self.data['TE']['raw_data'][0], self.data['TE']['raw_data'][1], color = cycle[0])
            self.ax5.plot(self.data['TE']['fitted_data'][0], self.data['TE']['fitted_data'][1], color = cycle[1])
            try:
                self.ax5.set_xlim((10, 0.9*1E9/self.data['TE']['metadata']['excitation_frequency']))
            except:
                self.ax5.set_xlim((10, 0.9*1E9/self.data['metadata']['TE']['excitation_frequency']))
            self.ax5.set_ylim((1E0,np.max(self.data['TE']['raw_data'][1])*2))

        self.canvas1.draw()


    # def update_table(self):
    #     self.tree.delete(*self.tree.get_children())
    #     peak_table = self.Chrom.create_peak_table(self.data['peak_data'], show_table=False, save_filename=None)
    #     for i, peak in enumerate(peak_table):
    #         if i != 0:
    #             self.tree.insert("", "end", values = (peak))

    def load_data(self, filename):

        with open(filename, 'rb') as f:
            self.data = pkl.load(f)
        #shoe filename
        self.label_2 = tk.Label(self, anchor = 'w', textvariable=self.text_filename)
        self.label_2.place(x=70, y = 10, height=20, width=400)
        self.text_filename.set('%s' %os.path.basename(filename))

        if 'job' in self.data.keys():
            sample_name = self.data['job']['target_name']
        elif 'sample' in self.data['metadata'].keys():
            sample_name = self.data['metadata']['sample']['name']
        else :
            sample_name = ''
        self.text_samplename.set(sample_name)
        if self.data['absorption'] is not None:
            self.text_absmax.set('{:.4f} ({:.1f} nm)'.format(self.data['absorption']['max'], self.data['absorption']['lambda_max']))
            self.text_absend.set('{:.1f} nm'.format(self.data['absorption']['end_wavelength']))
        if self.data['PL'] is not None:
            self.text_uvref.set('{:.4f} (mW)'.format(self.data['uv']['reference']*1000))
            self.text_uvabs.set('{:.4f} (mW)'.format(self.data['uv']['absorption']*1000))
            self.text_PLmax.set('{:.4f} ({:.1f} nm)'.format(self.data['PL']['max'], self.data['PL']['lambda_max']))
            self.text_relativeQY.set('{:.3f}'.format(self.data['PL']['relative_QY']))
            if self.data['PL']['max_gain_factor']:
                self.text_gainmax.set('{:.4f} ({:.1f} nm)'.format(self.data['PL']['max_gain_factor']*1E24, self.data['PL']['max_gain_wavelength']))
            if self.data['PL']['maintenance(at_1min)']:
                self.text_PLretention.set('{:.3f}'.format(self.data['PL']['maintenance(at_1min)']))
            if self.data['uv']['absorbance_maintenance(at_1min)']:
                self.text_Absretention.set('{:.3f}'.format(self.data['uv']['absorbance_maintenance(at_1min)']))
        if self.data['TE'] is not None:
            self.text_PLlifetime.set('{:.3f} ns'.format(self.data['TE']['fitting_results']['tau']))
        # self.update_table()
        self.update_plot()


    def select_data(self):

        data_num = self.r2_val.get()

        self.load_data(self.filenames[data_num])



    def file_open(self):
        self.filenames = filedialog.askopenfilenames(filetypes=[("pickle file", ".pkl")])

        self.load_data(self.filenames[0])

        if self.r_buttons:
            for i in range(len(self.r_buttons)):
                self.r_buttons[-i-1].destroy()

        self.r_buttons = []
        for i, filename in enumerate(self.filenames):
            fname = os.path.basename(filename)
            row = i%8
            column = i//8
            self.r_buttons.append(tk.Radiobutton(self, text = fname.split('.')[0], value = i, var = self.r2_val, command = self.select_data))
            self.r_buttons[i].place(x = 40 + 230 * column, y = 650 + 22 * row)


    def exit(self):
        self.master.quit()
        


 


    # def plot_mass_at_peak(self):
    #     self.get_peak_number()
    #     self.sub_win.wait_window()

    #     MS_spectra = self.data['DAD_data']['chromatograms']['peak_data'][self.peak_num]['MS_spectra']
    #     fig = plt.figure(figsize=(12,8))
    #     ax = []
    #     for i, (key, val) in enumerate(MS_spectra.items()):
    #         ax.append(fig.add_subplot(len(MS_spectra.keys()),1,i + 1, title = key))
    #         ax[i].bar(val[0],val[1])
    #     plt.show(block = False)


    # def plot_absorption_at_peak_all(self):

    #     self.Chrom.show_DAD_peak_spectra(self.data['DAD_data'], save_filename=None, block = False)

    # def plot_absorption_at_peak(self):
    #     self.get_peak_number()
    #     self.sub_win.wait_window()
    #     peak_data = self.data['DAD_data']['chromatograms']['peak_data'][self.peak_num]

    #     fig = plt.figure()
    #     ax = fig.add_subplot(1,1,1, title = 'Absorption spectrum at peak %s' %self.peak_num, xlabel = 'wavelength/nm', ylabel = 'absorbance/mAU')
    #     ax.plot(peak_data['spectrum'][0], peak_data['spectrum'][1])
    #     ax.legend(['{:.3f} min'.format(peak_data['peak']['time'])], loc = 'upper right')
    #     plt.show(block = False)

root = tk.Tk()
main = Application(master=root)
main.mainloop()


# def MS_analysis(filename, target_types = None, targets = None, target_dicts= None, target_plots = False, result_plot = True): 

#     if target_types:
#         target_list = []
#         for target_type in target_types:
#             target_list.extend(list(filter(lambda value:value['type'] == target_type , Chrom.target_dict.values())))

#     if targets:
#         target_list = [Chrom.target_dict[target] for target in targets]

#     if target_dicts:
#         target_list = target_dicts

#     #parameters for peak search for TIC and XIC
#     setting_TIC =   {'peak_deconvolution' : 'None', # None/gaussian/exponnorm/t
#                             'low_pass_cutoff' : 0.1, 
#                             'grad2_threshold' : 5, 
#                             'height_threshold_abs': 1E7,
#                             'height_threshold_rel' : 0.1,
#                             'grad_threshold_base' : 0.01, 
#                             'peak_split' : 1E6, 
#                             'bkg_corr' : True, 
#                             'ball_radius' : [1000, 20, 1E-8],
#                             'time_range' : [0.1, None], 
#                             'grad_plot' :  False,
#                             'bkg_plot' : False,
#                             'deconv_plot' : False,
#                             'result_plot' : False,
#                             'verbose'     : False 
#                             }  

#     setting_XIC =   {'peak_deconvolution' : 'None', # None/gaussian/exponnorm/t
#                             'low_pass_cutoff' : 0.05, 
#                             'grad2_threshold' : 20, 
#                             'height_threshold_abs': 3E5,
#                             'height_threshold_rel' : 0.05,
#                             'grad_threshold_base' : 0.03, 
#                             'peak_split' : 1E6, 
#                             'bkg_corr' : False, 
#                             'time_range' : [0.1, None], 
#                             'grad_plot' :  False,
#                             'bkg_plot' : False,
#                             'deconv_plot' : False,
#                             'result_plot' : False,
#                             'verbose'     : False
#                             }        

#     MSdata = Chrom.MS_analysis(filename, target_list, setting_TIC, setting_XIC,
#                         mass_tolerance= 0.002, istp_threshold = 0.5, target_plots = target_plots,  result_plot = result_plot) 

#     return MSdata

#print(injection.Name)
#Chrom.Device.show_MS_metadata(injection)
#Chrom.Device.add_injection('test2', sequence, 1, 'B:B1', instmeth, processing_method = None, s_type = 0, status = 0)
# #injection = Chrom.Device.load_injection_from_sequence(sequence_url, 2)
# injection = Chrom.load_injection(sequence_url, 2)
# print(injection.Name)
# print(Chrom.Device.get_gradient_table(injection.InstrumentMethod))
#Chrom.Device.change_script_time(instmeth, 'Run', 4, step_number = 2)

#Chrom.Device.list_instrument_symbol()
#Chrom.Device.show_PlugInData(inst_method)



###MS calculator###################################################

###create analysis data##
# DAD_filelist = filedialog.askopenfilenames()
# MS_filelist = filedialog.askopenfilenames()

# DAD_filelist = [1]

# for i, DAD_file in enumerate(DAD_filelist):     
# #     print(DAD_file, MS_filelist[i])
# #     DAD_result = Chrom.DAD_analysis_3D(DAD_file, **DAD_peak_search_setting)
# #     DAD_result = Chrom.add_MS_spectra_to_DAD_results(DAD_result, Chrom.load_MS_data(MS_filelist[i]), DAD_offset = -0.03, show_plot=False)

#      MS_data = MS_analysis(MS_filelist[i], 'BMIDA')

#     DAD_result = Chrom.DAD_peak_anotation_by_MSdata(DAD_result, MS_data, DAD_offset = -0.03, time_tolerance= 0.03, show_plot= False)

#     filename = 'D:/temp_for_covid19/HPLCMS/analysis_results/%s_res.pkl' %os.path.basename(DAD_file).split('.')[0]

#     Chrom.savefile(filename, DAD_result)

############

# result_file = 'D:/temp_for_covid19/HPLCMS/analysis_results/0_BA mix 10-100 f0_res.pkl'
# result_data = Chrom.loadfile(result_file)

# Chrom.create_DAD_peak_table(result_data, show_table=True)

# peak = result_data['chromatograms']['peak_data'][3]
# print(peak.keys())


# for key, val in peak['MS_spectra'].items():
#    plt.bar(val[0], val[1])
#    plt.title(key)
#    plt.xlim(80, 700)
#    plt.show()
# for i, peak in enumerate(result_data['chromatograms']['peak_data']):
#    print(i, peak['annotation'])


#filename = 'C:/Users/smily/Dropbox/PythonScript/kazu/data/formula_test.pkl'



# Chrom.savefile(filename, DAD_result)

# mass_spectrum = MSdata['Guaiazulene_XIC']['peak_data'][0]['spectrum']
#mass_spectrum = MSdata['6Methoxy2Naphtha_BA_XIC']['peak_data'][0]['spectrum']
#mass_spectrum = MSdata['Phenanthracene9_BA_XIC']['peak_data'][0]['spectrum']
#mass_spectrum = MSdata['2BrbiPhen2_BMIDA_XIC']['peak_data'][0]['spectrum']
#mass_spectrum = MSdata['5Br22biThiophen5_BMIDA_XIC']['peak_data'][0]['spectrum']
#mass_spectrum = MSdata['5CyThiophen2_BMIDA_XIC']['peak_data'][0]['spectrum']


# candidates = Chrom.mass_to_formula(399.9310874561875, 3, charge = -1, space = search_space, n_jobs = -2, N_rule='even')
# for candidate in candidates:
#    if 'C13' in candidate['mf']:
#        print(candidate)

#print(candidate_formula)
#print(len(candidate_formula))

#candidates = Chrom.guess_molecular_formula(mass_spectrum, 199.14812702009112, 5, charge = 1, space = search_space, istp_threshold = 0.5, n_jobs = -2, N_rule = 'even', do_plot = False)
#candidates = Chrom.guess_molecular_formula(mass_spectrum, 201.07284807990888, 3, charge = -1, space = search_space, istp_threshold = 0.5, n_jobs = -2, N_rule = 'even', do_plot = True)
#candidates = Chrom.guess_molecular_formula(mass_spectrum, 221.07793347990886, 3, charge = -1, space = search_space, istp_threshold = 0.5, n_jobs = -2, N_rule = 'even', do_plot = True)
#candidates = Chrom.guess_molecular_formula(mass_spectrum, 386.0204748799089, 3, charge = -1, space = search_space, istp_threshold = 0.5, n_jobs = -2, N_rule = 'even', do_plot = False)
#candidates = Chrom.guess_molecular_formula(mass_spectrum, 397.9328 + 5.4857990888E-4, 3, charge = -1, space = search_space, istp_threshold = 0.3, n_jobs = -2, N_rule = 'even', do_plot = False)

#for candidate in candidates:
#   print(candidate)

###########

# search_space = 'C2-100 H0-200 N0-10 O0-10 B0-3 F0-6 Cl0-6 Br0-6 I0-6 S0-3'

# filename = 'D:/temp_for_covid19/HPLCMS/analysis_results/0_MIDA mix 10-100 f0_res.pkl'
# DAD_result = Chrom.loadfile(filename)

# # Chrom.create_DAD_peak_table(DAD_result, show_table=True)

# candidate = Chrom.guess_molecular_formula_at_DAD_peak(DAD_result['chromatograms']['peak_data'][20], tolerance = 5, space = search_space, istp_threshold = 0.5, n_jobs = -2, N_rule = 'even', do_plot = False)
# for key, val in candidate.items():
#     print(key)
#     for v in val:
#         print(v)
############

#BMIDA: C5H7O4NB
#C13H11S2BrO4NB
#C8H4S2BrBmida