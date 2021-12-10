from matplotlib.figure import Figure
from tkinter import *
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import pprint as p
import ECSim as ecs
from functools import partial
import pickle
from tkinter import filedialog
import pandas as pd


FONT = ('Arial', 12)


class mclass:

    def __init__(self, window):

        ################################################################################################################
        # General Setup ################################################################################################
        ################################################################################################################

        self.window = window

        width = window.winfo_screenwidth()
        height = window.winfo_screenheight()
        window.geometry("%dx%d" % (width, height))

        self.Var_Frame = LabelFrame(window, text='Variables', font=FONT, padx=5, pady=5)
        self.Var_Frame.grid(row=0, column=0, padx=10, pady=10)

        self.menubar = Menu(self.window,
                            background='#ff8000', foreground='black', activebackground='white',
                            activeforeground='black')
        self.file = Menu(self.menubar, tearoff=1, background='#ffcc99', foreground='black')
        self.file.add_command(label="New")
        self.file.add_command(label="Open", command=self.load_data)
        self.file.add_command(label="Save as", command=self.save_data)
        self.file.add_separator()
        self.file.add_command(label="Exit", command=window.quit)
        self.menubar.add_cascade(label="File", menu=self.file)
        self.menubar.config(font=FONT)
        self.window.config(menu=self.menubar)

        self.DataDic = {}

        self.MultiPlot = {
            'Numplot': 3,
            'varType': 'SR',
            'varType2': '',
            'var_Low': 0.01,
            'var_High': 0.0001,
            'Label': '',
            'varValue': 0,
            'varTest': 0,
            'Plots': {},
            'Active': False,
        }

        self.Parameters = {
            'init_E': -0.5,
            'high_E': 1.7,
            'low_E': -0.5,
            'final_E': -0.5,
            'enableFinal': False,
            'direction': 'pos',
            'SR': 0.1,
            'seg': 2,
            'segList': [],
            'eType': 'disk',
            'eSize': 3.0e-3,
        }

        self.Solution = {}
        self.Electron_Step = {}
        self.Chemical_Step = {}

        # Save / Load
        self.rxnDic = []
        self.rdxDic = []

        # Buttons
        self.button = Button(window, text="Plot Surface", font=FONT, command=self.plot_IBL)
        self.button.grid(row=3, column=0, sticky="n", padx=5, pady=5)

        self.button2 = Button(window, text="Export Plot", font=FONT, command=self.export_data)
        self.button2.grid(row=4, column=0, sticky="n", padx=5, pady=5)

        ################################################################################################################
        # MultiPlot UI #################################################################################################
        ################################################################################################################

        self.MultiPlot_Active = True
        self.multiDic = {}

        # Frame
        self.MultiPlot_Frame = LabelFrame(self.Var_Frame, font=FONT, text='Multiple Plots', padx=5, pady=5)
        self.MultiPlot_Frame.grid(row=0, column=0, padx=10, pady=10)

        # Dropdown
        Numplot_Options = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10']
        self.Drop_Numplot_Var = StringVar()
        self.Drop_Numplot_Var.set(Numplot_Options[0])
        self.Drop_Numplot = OptionMenu(self.MultiPlot_Frame, self.Drop_Numplot_Var, *Numplot_Options)
        self.Drop_Numplot.grid(row=22, column=1, padx=2, pady=2)
        self.L_Drop_Numplot = Label(self.MultiPlot_Frame, font=FONT, text='Number of Plots:')
        self.L_Drop_Numplot.grid(row=22, column=0, padx=2, pady=2)

        varType_Options = ['SR', 'k_f', 'k_b', 'alpha']
        self.Drop_varType_Var = StringVar()
        self.Drop_varType_Var.set(varType_Options[0])
        self.Drop_varType = OptionMenu(self.MultiPlot_Frame, self.Drop_varType_Var, *varType_Options,
                                       command=self.refresh_varType)
        self.Drop_varType.grid(row=23, column=1, padx=2, pady=2)
        self.L_Drop_varType = Label(self.MultiPlot_Frame, font=FONT, text='Modifying:')
        self.L_Drop_varType.grid(row=23, column=0, padx=2, pady=2)

        varType2_Options = ['None']
        self.Drop_varType2_Var = StringVar()
        self.Drop_varType2_Var.set(varType2_Options[0])
        self.Drop_varType2 = OptionMenu(self.MultiPlot_Frame, self.Drop_varType2_Var, *varType2_Options)
        self.Drop_varType2.grid(row=24, column=1, padx=2, pady=2)
        self.L_Drop_varType2 = Label(self.MultiPlot_Frame, font=FONT, text='rdx/rxn:')
        self.L_Drop_varType2.grid(row=24, column=0, padx=2, pady=2)

        # Input Fields
        self.multiDic['Lowest Value_var'] = DoubleVar()
        self.multiDic['Highest Value_var'] = DoubleVar()
        DD = {}
        n = 25
        for keys in self.multiDic.keys():
            i = keys.replace('_var', '')
            DD[i] = Entry(self.MultiPlot_Frame, width=10, font=FONT, textvariable=self.multiDic[keys]).grid(
                row=n, column=1, padx=2, pady=2)
            DD['L_' + i] = Label(self.MultiPlot_Frame, font=FONT, text=i + ':').grid(row=n, column=0, padx=2, pady=2)
            n += 1
        temp = {**DD, **self.multiDic}
        self.multiDic = temp
        del DD
        del temp
        del n

        # Buttons/Checkboxes
        self.mp_show = Button(self.Var_Frame, font=FONT, text='Activate MultiPlot', command=self.showMP)
        self.mp_show.grid(row=1, column=0, padx=2, pady=2)
        self.mp_hide = Button(self.Var_Frame, font=FONT, text='Disable MultiPlot', command=self.hideMP)
        self.mp_hide.grid(row=1, column=0, padx=2, pady=2)
        self.mp_show.grid_remove()

        ################################################################################################################
        # CV Parameters UI #############################################################################################
        ################################################################################################################

        self.paramDic = {}

        # Frame
        self.Parameters_Frame = LabelFrame(self.Var_Frame, font=FONT, text='CV Parameters', padx=5, pady=5)
        self.Parameters_Frame.grid(row=2, column=0, padx=10, pady=10)

        # Dropdown
        direction_Options = ['pos', 'neg']
        self.Drop_direction_Var = StringVar()
        self.Drop_direction_Var.set(direction_Options[0])
        self.Drop_direction = OptionMenu(self.Parameters_Frame, self.Drop_direction_Var, *direction_Options)
        self.Drop_direction.grid(row=24, column=1, padx=2, pady=2)
        self.L_Drop_direction = Label(self.Parameters_Frame, font=FONT, text='Direction:')
        self.L_Drop_direction.grid(row=24, column=0, padx=2, pady=2)

        eType_Options = ['disk', 'cylinder', 'hemisphere', 'rect', 'square', 'sphere']
        self.Drop_eType_Var = StringVar()
        self.Drop_eType_Var.set(eType_Options[0])
        self.Drop_eType = OptionMenu(self.Parameters_Frame, self.Drop_eType_Var, *eType_Options)
        self.Drop_eType.grid(row=69, column=1, padx=2, pady=2)
        self.L_Drop_eType = Label(self.Parameters_Frame, font=FONT, text='Electrode Type:')
        self.L_Drop_eType.grid(row=69, column=0, padx=2, pady=2)

        # Input Fields
        self.paramDic['Initial E_var'] = DoubleVar()
        self.paramDic['High E_var'] = DoubleVar()
        self.paramDic['Low E_var'] = DoubleVar()
        self.paramDic['Final E_var'] = DoubleVar()
        self.paramDic['Scan Rate (V/s)_var'] = DoubleVar()
        self.paramDic['Segments_var'] = IntVar()
        self.paramDic['Electrode Size (m)_var'] = DoubleVar()
        DD = {}
        n = 25
        for keys in self.paramDic.keys():
            i = keys.replace('_var', '')
            DD[i] = Entry(self.Parameters_Frame, width=10, font=FONT, textvariable=self.paramDic[keys]).grid(
                row=n, column=1, padx=2, pady=2)
            DD['L_' + i] = Label(self.Parameters_Frame, font=FONT, text=i + ':').grid(row=n, column=0, padx=2, pady=2)
            n += 1
        temp = {**DD, **self.paramDic}
        self.paramDic = temp
        del DD
        del temp
        del n

        # Buttons/Checkboxes
        self.Check_enableFinal_Var = BooleanVar()
        self.Check_enableFinal = Checkbutton(self.Parameters_Frame,
                                             variable=self.Check_enableFinal_Var, onvalue=True, offvalue=False)
        self.Check_enableFinal.grid(row=70, column=1, padx=2, pady=2)
        self.L_Check_enableFinal = Label(self.Parameters_Frame, font=FONT, text='Enable Final E:')
        self.L_Check_enableFinal.grid(row=70, column=0, padx=2, pady=2)

        ################################################################################################################
        # Redox Parameters UI ##########################################################################################
        ################################################################################################################

        self.rdxList = []

        # Adding new rdx conditions
        self.rdx_canvas = Canvas(window, borderwidth=0, width=250, height=500, highlightthickness=1)
        self.rdx_lframe = LabelFrame(self.rdx_canvas, font=FONT, text="Redox Parameters")
        self.rdx_sbar = Scrollbar(window, orient="vertical", command=self.rdx_canvas.yview)
        self.rdx_canvas.configure(yscrollcommand=self.rdx_sbar.set)

        self.rdx_sbar.grid(row=0, column=1, rowspan=2, sticky="nsw")
        self.rdx_canvas.grid(row=0, column=2, rowspan=2, sticky="nsew")
        self.rdx_canvas.create_window((2, 2), window=self.rdx_lframe, anchor="n",
                                      tags="self.rdx_lframe")

        self.rdx_lframe.bind("<Configure>", self.rdx_onFrameConfigure)
        self.rdx_rframe = Frame(window)
        self.rdx_rframe.grid(row=2, column=2, sticky="n")

        rdx_add_btn = Button(self.rdx_rframe, font=FONT, text='Add Redox Event', command=self.add_rdx)
        rdx_add_btn.grid(row=1, column=0, sticky="nwe")
        rdx_remove_btn = Button(self.rdx_rframe, font=FONT, text='Remove Redox Event', command=self.remove_rdx)
        rdx_remove_btn.grid(row=2, column=0, sticky="nwe")

        self.rdx_cnt = 1

        ################################################################################################################
        # Reaction Parameters UI #######################################################################################
        ################################################################################################################

        self.rxnList = []
        # Adding new rxn conditions
        self.rxn_canvas = Canvas(window, borderwidth=0, width=380, height=400, highlightthickness=0)
        self.rxn_lframe = LabelFrame(self.rxn_canvas, font=FONT, text="Reaction Parameters")
        self.rxn_sbar = Scrollbar(window, orient="vertical", command=self.rxn_canvas.yview)
        self.rxn_canvas.configure(yscrollcommand=self.rxn_sbar.set)

        self.rxn_sbar.grid(row=0, column=3, rowspan=2, sticky="nsw")
        self.rxn_canvas.grid(row=0, column=4, rowspan=2, sticky="nsew")
        self.rxn_canvas.create_window((2, 2), window=self.rxn_lframe, anchor="n",
                                      tags="self.rxn_lframe")

        self.rxn_lframe.bind("<Configure>", self.rxn_onFrameConfigure)
        self.rxn_rframe = Frame(window)
        self.rxn_rframe.grid(row=2, column=4, sticky="nw")

        rxn_add_btn = Button(self.rxn_rframe, font=FONT, text='Add Reaction Event', command=self.add_rxn)
        rxn_add_btn.grid(row=1, column=0, sticky="nwe")
        rxn_remove_btn = Button(self.rxn_rframe, font=FONT, text='Remove Reaction Event', command=self.remove_rxn)
        rxn_remove_btn.grid(row=2, column=0, sticky="nwe")
        rxn_refresh_btn = Button(self.rxn_rframe, font=FONT, text='Refresh Reaction Event',
                                 command=self.refresh_rxn_options)
        rxn_refresh_btn.grid(row=1, column=1, sticky="nwe")
        # getList_btn = Button(self.rxn_rframe, font=FONT, text='Print Lists', command=self.getLists)
        # getList_btn.grid(row=2, column=1, sticky='nwe')
        # Save_btn = Button(self.rxn_rframe, font=FONT, text='Save', command=self.save_data)
        # Save_btn.grid(row=3, column=0, sticky='nwe')
        # Load_btn = Button(self.rxn_rframe, font=FONT, text='Load', command=self.load_data)
        # Load_btn.grid(row=3, column=1, sticky='nwe')

        self.rxn_cnt = 1

        ################################################################################################################
        # Graphing UI ##################################################################################################
        ################################################################################################################

        self.graph_canvas = None
        self.Plot_Frame = LabelFrame(window, padx=5, pady=5)
        self.Plot_Frame.grid(row=0, column=1, padx=10, pady=10)

    ########################################################################################################################

    def add_rdx(self):

        rdx_frame = LabelFrame(self.rdx_lframe, font=FONT, text="rdx" + str(self.rdx_cnt))
        rdx_frame.grid(row=self.rdx_cnt, column=0, padx=5, pady=5)

        self.rdxList.append({})
        n = len(self.rdxList) - 1
        # add variables to rdxList
        self.rdxList[n]['nameOx_var'] = StringVar()
        self.rdxList[n]['concOx_var'] = DoubleVar()
        self.rdxList[n]['diffOx_var'] = DoubleVar()
        self.rdxList[n]['nameRed_var'] = StringVar()
        self.rdxList[n]['concRed_var'] = DoubleVar()
        self.rdxList[n]['diffRed_var'] = DoubleVar()
        self.rdxList[n]['n_var'] = IntVar()
        self.rdxList[n]['E_0_var'] = DoubleVar()
        self.rdxList[n]['k_e_var'] = DoubleVar()
        self.rdxList[n]['alpha_var'] = DoubleVar()

        k = 0
        DD = {}
        temp = {}
        for i in self.rdxList[n].keys():
            if '_var' in i:
                j = i.replace('_var', '')
                DD[j] = Entry(rdx_frame, width=15, font=FONT, textvariable=self.rdxList[n][i])
                DD['L_' + j] = Label(rdx_frame, font=FONT, text=[j + ':'])
                DD[j].grid(row=k, column=1, padx=2, pady=2)
                DD['L_' + j].grid(row=k, column=0, padx=2, pady=2, sticky='e')
                k += 1
        temp = {**DD, **self.rdxList[n]}
        self.rdxList[n] = temp
        del DD
        del temp
        # p.pprint(self.rdxList)
        # print("\n")
        # print(self.rdxList[n]['nameOx'].master() is rdx_frame)
        self.rdx_cnt += 1

    ########################################################################################################################

    def remove_rdx(self):
        if len(self.rdxList) != 0:
            tempMaster = self.rdxList[-1]['nameOx'].master
            for i in self.rdxList[-1].keys():
                if '_var' not in i:
                    for widget in self.rdxList[-1][i].master.winfo_children():
                        widget.destroy()
            tempMaster.destroy()
            del self.rdxList[-1]
            self.rdx_cnt -= 1
            # p.pprint(self.rdxList)
            # print("\n")

    ########################################################################################################################

    def add_rxn(self):

        rxn_frame = LabelFrame(self.rxn_lframe, font=FONT, text="rxn" + str(self.rxn_cnt))
        rxn_frame.grid(row=self.rxn_cnt, column=0, padx=5, pady=5)
        rxn_frame2 = Frame(rxn_frame)
        rxn_frame2.grid(row=0, column=0, padx=5, pady=5, columnspan=2)

        self.rxnList.append({})
        n = len(self.rxnList) - 1
        # add variables to rxnList
        self.rxnList[n]['nameA_var'] = StringVar()
        self.rxnList[n]['concA_var'] = DoubleVar()
        self.rxnList[n]['diffA_var'] = DoubleVar()

        self.rxnList[n]['nameB_var'] = StringVar()
        self.rxnList[n]['concB_var'] = DoubleVar()
        self.rxnList[n]['diffB_var'] = DoubleVar()

        self.rxnList[n]['nameC_var'] = StringVar()
        self.rxnList[n]['concC_var'] = DoubleVar()
        self.rxnList[n]['diffC_var'] = DoubleVar()

        self.rxnList[n]['nameD_var'] = StringVar()
        self.rxnList[n]['concD_var'] = DoubleVar()
        self.rxnList[n]['diffD_var'] = DoubleVar()
        self.rxnList[n]['k_f_var'] = DoubleVar()
        self.rxnList[n]['k_b_var'] = DoubleVar()

        k = 5
        DD = {}
        temp = {}
        for i in self.rxnList[n].keys():
            if '_var' in i:
                j = i.replace('_var', '')
                DD[j] = Entry(rxn_frame, width=20, font=FONT, textvariable=self.rxnList[n][i])
                DD['L_' + j] = Label(rxn_frame, font=FONT, text=[j])
                DD[j].grid(row=k, column=1, padx=2, pady=2)
                DD['L_' + j].grid(row=k, column=0, padx=2, pady=2, sticky='e')
                k += 1
        temp = {**DD, **self.rxnList[n]}
        self.rxnList[n] = temp
        del DD
        del temp

        self.rxnList[n]['Options'] = ['None']
        """
        for rdxnum in range(len(self.rdxList)):
            for name in self.rdxList[rdxnum].keys():
                if 'name' in name:
                    if 'L_' not in name:
                        if '_var' not in name:
                            self.rxnList[n]['Options'].append(self.rdxList[rdxnum][name+'_var'].get())
        """
        # print('\n')
        # print(self.rxnList[n]['Options'])
        # print('\n')

        self.rxnList[n]['Drop1_var'] = StringVar()
        self.rxnList[n]['Drop1_var'].set(self.rxnList[n]['Options'][0])
        self.rxnList[n]['Drop2_var'] = StringVar()
        self.rxnList[n]['Drop2_var'].set(self.rxnList[n]['Options'][0])
        self.rxnList[n]['Drop3_var'] = StringVar()
        self.rxnList[n]['Drop3_var'].set(self.rxnList[n]['Options'][0])
        self.rxnList[n]['Drop4_var'] = StringVar()
        self.rxnList[n]['Drop4_var'].set(self.rxnList[n]['Options'][0])
        self.rxnList[n]['Drop1'] = OptionMenu(rxn_frame2, self.rxnList[n]['Drop1_var'], *self.rxnList[n]['Options'])
        self.rxnList[n]['Drop2'] = OptionMenu(rxn_frame2, self.rxnList[n]['Drop2_var'], *self.rxnList[n]['Options'])
        self.rxnList[n]['Drop3'] = OptionMenu(rxn_frame2, self.rxnList[n]['Drop3_var'], *self.rxnList[n]['Options'])
        self.rxnList[n]['Drop4'] = OptionMenu(rxn_frame2, self.rxnList[n]['Drop4_var'], *self.rxnList[n]['Options'])
        self.rxnList[n]['Drop1'].grid(row=0, column=0, padx=2, pady=2)
        self.rxnList[n]['Drop2'].grid(row=0, column=2, padx=2, pady=2)
        self.rxnList[n]['Drop3'].grid(row=2, column=0, padx=2, pady=2)
        self.rxnList[n]['Drop4'].grid(row=2, column=2, padx=2, pady=2)
        plus1 = Label(rxn_frame2, font=FONT, text='+').grid(row=0, column=1, padx=2, pady=2)
        plus2 = Label(rxn_frame2, font=FONT, text='+').grid(row=2, column=1, padx=2, pady=2)
        arrow = Label(rxn_frame2, font=FONT, text='-->').grid(row=1, column=1, padx=2, pady=2)

        # p.pprint(self.rxnList)
        # print("\n")
        # print(self.rxnList[n]['nameOx'].master() is rxn_frame)
        self.rxn_cnt += 1

    ########################################################################################################################

    def remove_rxn(self):
        if len(self.rxnList) != 0:
            tempMaster = self.rxnList[-1]['nameA'].master
            """for i in self.rxnList[-1].keys():
                if '_var' not in i:
                    if 'Options' not in i:
                        for widget in self.rxnList[-1][i].master.winfo_children():
                            widget.destroy()"""
            tempMaster.destroy()
            del self.rxnList[-1]
            self.rxn_cnt -= 1
            # p.pprint(self.rxnList)
            # print("\n")

    ########################################################################################################################

    def refresh_rxn_options(self):
        for rxnnum in range(len(self.rxnList)):
            self.rxnList[rxnnum]['Options'] = ['None']
            for rdxnum in range(len(self.rdxList)):
                for name in self.rdxList[rdxnum].keys():
                    if 'name' in name:
                        if 'L_' not in name:
                            if '_var' not in name:
                                if self.rdxList[rdxnum][name + '_var'].get() not in '':
                                    self.rxnList[rxnnum]['Options'].append(self.rdxList[rdxnum][name + '_var'].get())

            for name in self.rxnList[rxnnum].keys():
                if 'name' in name:
                    if 'L_' not in name:
                        if '_var' not in name:
                            if self.rxnList[rxnnum][name + '_var'].get() not in '':
                                self.rxnList[rxnnum]['Options'].append(self.rxnList[rxnnum][name + '_var'].get())
            for DropNum in range(4):
                # print(DropNum+1, 'DropNumber!')
                # print('\n')
                menu = self.rxnList[rxnnum]['Drop' + str(DropNum + 1)]['menu']
                menu.delete(0, 'end')
                for string in self.rxnList[rxnnum]['Options']:
                    menu.add_command(label=string,
                                     command=partial(self.rxnList[rxnnum][
                                                         'Drop' + str(DropNum + 1) + '_var'].set, string))
                    # print(string, 'Drop' + str(DropNum + 1) + '_var')
        # print('\n')
        # print(self.rxnList)
        # print('\n')

    ########################################################################################################################

    def refresh_varType(self, var):

        rdx_temp = []
        rxn_temp = []
        for rdx in range(len(self.rdxList)):
            rdx_temp.append('rdx' + str(rdx + 1))
        for rxn in range(len(self.rxnList)):
            rxn_temp.append('rxn' + str(rxn + 1))

        varType2_Options = ['None']
        varType2_Options.extend(rdx_temp)
        varType2_Options.extend(rxn_temp)

        menu = self.Drop_varType2['menu']
        menu.delete(0, 'end')
        for string in varType2_Options:
            menu.add_command(label=string,
                             command=partial(self.Drop_varType2_Var.set, string))

    ########################################################################################################################

    def rdx_onFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.rdx_canvas.configure(scrollregion=self.rdx_canvas.bbox("all"))

    ########################################################################################################################

    def rxn_onFrameConfigure(self, event):
        '''Reset the scroll region to encompass the inner frame'''
        self.rxn_canvas.configure(scrollregion=self.rxn_canvas.bbox("all"))

    ########################################################################################################################

    def showMP(self):
        self.MultiPlot_Frame.grid()
        self.mp_show.grid_remove()
        self.mp_hide.grid()
        self.MultiPlot_Active = True

    ########################################################################################################################

    def hideMP(self):
        self.MultiPlot_Frame.grid_remove()
        self.mp_hide.grid_remove()
        self.mp_show.grid()
        self.MultiPlot_Active = False

    ########################################################################################################################

    def clear_space(self):  # new
        self.graph_canvas._tkcanvas.destroy()
        # global main, root
        # main.destroy()
        # main = Frame(root)
        # main.pack()

    ########################################################################################################################

    def plot_IBL(self):

        self.DataDic = {}

        self.MultiPlot = {
            'Numplot': 3,
            'varType': 'SR',
            'var_Low': 0.01,
            'var_High': 0.0001,
            'Label': '',
            'varValue': 0,
            'varTest': 0,
            'Plots': {},
            'Active': False,
        }

        self.Parameters = {
            'init_E': -0.5,
            'high_E': 1.7,
            'low_E': -0.5,
            'final_E': -0.5,
            'enableFinal': False,
            'direction': 'pos',
            'SR': 0.1,
            'seg': 2,
            'segList': [],
            'eType': 'disk',
            'eSize': 3.0e-3,
        }

        if self.graph_canvas is not None:
            self.clear_space()
        Graph_Frame = LabelFrame(self.window, padx=5, pady=5)
        Graph_Frame.grid(row=0, column=50, padx=10, pady=10)
        fig = Figure(figsize=(8, 5))
        self.graph_canvas = FigureCanvasTkAgg(fig, master=Graph_Frame)
        self.graph_canvas.get_tk_widget().pack()
        self.graph_canvas._tkcanvas.pack(side="top", fill="both", expand=1)
        a = fig.add_subplot(111)
        a.set_xlabel("Potential (V)")
        a.set_ylabel("Current (A)")

        ######################
        # Pull Data and Plot #
        ######################

        # MultiPlot Dictionary
        self.MultiPlot['Numplot'] = int(self.Drop_Numplot_Var.get())
        self.MultiPlot['varType'] = self.Drop_varType_Var.get()
        self.MultiPlot['varType2'] = self.Drop_varType2_Var.get()
        self.MultiPlot['var_Low'] = self.multiDic['Lowest Value_var'].get()
        self.MultiPlot['var_High'] = self.multiDic['Highest Value_var'].get()
        self.MultiPlot['Active'] = self.MultiPlot_Active

        print('\n')
        print('MultiPlot')
        p.pprint(self.MultiPlot)

        # Parameters Dictionary

        self.Parameters['init_E'] = self.paramDic['Initial E_var'].get()
        self.Parameters['high_E'] = self.paramDic['High E_var'].get()
        self.Parameters['low_E'] = self.paramDic['Low E_var'].get()
        self.Parameters['final_E'] = self.paramDic['Final E_var'].get()
        self.Parameters['enableFinal'] = self.Check_enableFinal_Var.get()
        self.Parameters['direction'] = self.Drop_direction_Var.get()
        self.Parameters['SR'] = self.paramDic['Scan Rate (V/s)_var'].get()
        self.Parameters['seg'] = self.paramDic['Segments_var'].get()
        self.Parameters['eType'] = self.Drop_eType_Var.get()
        self.Parameters['eSize'] = self.paramDic['Electrode Size (m)_var'].get()

        print('\n')
        print('CV Parameters')
        p.pprint(self.Parameters)

        # Solution Dictionary (slightly different format)
        for rdxNum in range(len(self.rdxList)):
            for x in range(2):
                if x == 0:
                    if self.rdxList[rdxNum]['nameOx_var'].get() != '':
                        self.Solution['Ox' + str(rdxNum)] = [
                            ecs.ecs.Species(self.rdxList[rdxNum]['nameOx_var'].get(),
                                            self.rdxList[rdxNum]['concOx_var'].get(),
                                            self.rdxList[rdxNum]['diffOx_var'].get(),
                                            ),
                            self.rdxList[rdxNum]['nameOx_var'].get(),
                        ]
                else:
                    if self.rdxList[rdxNum]['nameOx_var'].get() != '':
                        self.Solution['Red' + str(rdxNum)] = [
                            ecs.ecs.Species(self.rdxList[rdxNum]['nameRed_var'].get(),
                                            self.rdxList[rdxNum]['concRed_var'].get(),
                                            self.rdxList[rdxNum]['diffRed_var'].get(),
                                            ),
                            self.rdxList[rdxNum]['nameRed_var'].get(),
                        ]
        for rxnNum in range(len(self.rxnList)):
            for keys in self.rxnList[rxnNum]:
                for x in ['A', 'B', 'C', 'D']:
                    if self.rxnList[rxnNum]['name' + x].get() != '':
                        self.Solution[x + str(rxnNum)] = [
                            ecs.ecs.Species(self.rxnList[rxnNum]['name' + x + '_var'].get(),
                                            self.rxnList[rxnNum]['conc' + x + '_var'].get(),
                                            self.rxnList[rxnNum]['diff' + x + '_var'].get(),
                                            ),
                            self.rxnList[rxnNum]['name' + x + '_var'].get(),
                        ]

        print('\n')
        print('Solution')
        p.pprint(self.Solution)

        # Electron_Step Dictionary

        for rdxNum in range(len(self.rdxList)):
            if self.rdxList[rdxNum]['nameOx_var'].get() != '':
                self.Electron_Step['rdx' + str(rdxNum + 1)] = ecs.ecs.Redox(self.Solution['Ox' + str(rdxNum)][0],
                                                                            self.Solution['Red' + str(rdxNum)][0],
                                                                            self.rdxList[rdxNum]['n_var'].get(),
                                                                            self.rdxList[rdxNum]['E_0_var'].get(),
                                                                            self.rdxList[rdxNum]['k_e_var'].get(),
                                                                            self.rdxList[rdxNum]['alpha_var'].get(),
                                                                            ).enable()

        print('\n')
        print('Electron Step')
        p.pprint(self.Electron_Step)

        # Chemical Step (hard one)

        for rxnNum in range(len(self.rxnList)):
            A = self.rxnList[rxnNum]['Drop1_var'].get()
            B = self.rxnList[rxnNum]['Drop2_var'].get()
            C = self.rxnList[rxnNum]['Drop3_var'].get()
            D = self.rxnList[rxnNum]['Drop4_var'].get()
            vars = [A, B, C, D]
            for var in range(len(vars)):
                if vars[var] == 'None':
                    vars[var] = None
                else:
                    for keys in self.Solution.keys():
                        if vars[var] == self.Solution[keys][1]:
                            vars[var] = self.Solution[keys][0]

            self.Chemical_Step['rxn' + str(rxnNum + 1)] = ecs.ecs.Reaction(vars[0],
                                                                           vars[1],
                                                                           vars[2],
                                                                           vars[3],
                                                                           self.rxnList[rxnNum]['k_f_var'].get(),
                                                                           self.rxnList[rxnNum]['k_b_var'].get()
                                                                           ).enable()

        print('\n')
        print('Chemical Step')
        p.pprint(self.Chemical_Step)
        print('\n')

        y = ecs.VarSend(self.MultiPlot, self.Parameters, self.Electron_Step, self.Chemical_Step)
        self.DataDic = y[1]
        # p.pprint(self.DataDic)
        y = y[0]
        for plots in y.keys():
            a.plot(y[plots][0], y[plots][1], label=y[plots][2])
        a.legend(loc="upper left")
        a.grid('on', linestyle='--')
        self.graph_canvas.draw()
        self.graph_canvas.get_tk_widget().pack()
        toolbar = NavigationToolbar2Tk(self.graph_canvas, Graph_Frame)
        toolbar.update()
        self.graph_canvas.get_tk_widget().pack()
        width = self.window.winfo_screenwidth()
        height = self.window.winfo_screenheight()
        self.window.geometry("%dx%d" % (width, height))

        self.MultiPlot = {}
        self.Parameters = {}
        self.Solution = {}

    ########################################################################################################################

    def addField(self):
        self.rdxList.append({})
        n = len(self.rdxList) - 1

    ########################################################################################################################

    def getLists(self):

        print('\n')
        print('Rdx List:')
        p.pprint(self.rdxList)
        print('\n')
        print('Rxn List:')
        p.pprint(self.rxnList)
        print('\n')

    ########################################################################################################################

    def save_data(self, file_path=None):
        self.MultiPlot = {}
        self.Parameters = {}
        self.rdxDic = []
        self.rxnDic = []
        try:
            # rdx save
            for n in range(len(self.rdxList)):
                self.rdxDic.append({})
                for keys in self.rdxList[n].keys():
                    if '_var' in keys:
                        self.rdxDic[n][keys] = self.rdxList[n][keys].get()
            # rxn save
            for n in range(len(self.rxnList)):
                self.rxnDic.append({})
                self.rxnDic[n]['Options'] = self.rxnList[n]['Options']
                for keys in self.rxnList[n].keys():
                    if '_var' in keys:
                        self.rxnDic[n][keys] = self.rxnList[n][keys].get()
            # multi save
            if self.MultiPlot_Active:
                self.MultiPlot['Numplot'] = int(self.Drop_Numplot_Var.get())
                self.MultiPlot['varType'] = self.Drop_varType_Var.get()
                self.MultiPlot['varType2'] = self.Drop_varType2_Var.get()
                self.MultiPlot['var_Low'] = self.multiDic['Lowest Value_var'].get()
                self.MultiPlot['var_High'] = self.multiDic['Highest Value_var'].get()
                self.MultiPlot['Active'] = self.MultiPlot_Active

            # param save
            self.Parameters['init_E'] = self.paramDic['Initial E_var'].get()
            self.Parameters['high_E'] = self.paramDic['High E_var'].get()
            self.Parameters['low_E'] = self.paramDic['Low E_var'].get()
            self.Parameters['final_E'] = self.paramDic['Final E_var'].get()
            self.Parameters['enableFinal'] = self.Check_enableFinal_Var.get()
            self.Parameters['direction'] = self.Drop_direction_Var.get()
            self.Parameters['SR'] = self.paramDic['Scan Rate (V/s)_var'].get()
            self.Parameters['seg'] = self.paramDic['Segments_var'].get()
            self.Parameters['eType'] = self.Drop_eType_Var.get()
            self.Parameters['eSize'] = self.paramDic['Electrode Size (m)_var'].get()

            data = {
                'multi': self.MultiPlot,
                'param': self.Parameters,
                'rdx': self.rdxDic,
                'rxn': self.rxnDic,
            }
            if file_path is None:
                file_path = filedialog.asksaveasfilename(
                    defaultextension='.pickle', filetypes=[('Pickle files', '*.pickle'), ('All files', '*')]
                )

            with open(file_path, 'wb') as f:
                pickle.dump(data, f)
            print('Save Done')

        except Exception as e:
            print('error saving state:', str(e))

    ########################################################################################################################

    def export_data(self, file_path=None):
        if file_path is None:
            file_path = filedialog.asksaveasfilename(defaultextension='.csv', filetypes=[('csv files', '*.csv')])
        df = pd.DataFrame(self.DataDic)
        df.to_csv(file_path, index=False, header=True)
        print('Export Done')

    ########################################################################################################################

    def load_data(self, file_path=None):

        self.MultiPlot = {}
        self.Parameters = {}
        self.rdxDic = []
        self.rxnDic = []
        try:
            if file_path is None:
                file_path = filedialog.askopenfilename(
                    defaultextension='.pickle', filetypes=[('Pickle files', '*.pickle'), ('All files', '*')]
                )

            with open(file_path, 'rb') as f:
                data = pickle.load(f)

            self.MultiPlot = data['multi']
            self.Parameters = data['param']
            self.rdxDic = data['rdx']
            self.rxnDic = data['rxn']

            # rdx Load
            if len(self.rdxDic) == len(self.rdxList):
                pass
            else:
                diff = len(self.rdxDic) - len(self.rdxList)
                if diff > 0:
                    for n in range(diff):
                        self.add_rdx()
                        print('add:', n)
                else:
                    for n in range(-diff):
                        self.remove_rdx()
                        print('remove:', n)

            for n in range(len(self.rdxDic)):
                for keys in self.rdxDic[n].keys():
                    self.rdxList[n][keys].set(self.rdxDic[n][keys])

            # rxn Load
            if len(self.rxnDic) == len(self.rxnList):
                pass
            else:
                diff = len(self.rxnDic) - len(self.rxnList)
                if diff > 0:
                    for n in range(diff):
                        self.add_rxn()
                        print('add:', n)
                else:
                    for n in range(diff):
                        self.remove_rxn()
                        print('remove:', n)

            for n in range(len(self.rxnDic)):
                for keys in self.rxnDic[n].keys():
                    if 'Options' not in keys:
                        self.rxnList[n][keys].set(self.rxnDic[n][keys])

            # Param Load
            self.paramDic['Initial E_var'].set(self.Parameters['init_E'])
            self.paramDic['High E_var'].set(self.Parameters['high_E'])
            self.paramDic['Low E_var'].set(self.Parameters['low_E'])
            self.paramDic['Final E_var'].set(self.Parameters['final_E'])
            self.Check_enableFinal_Var.set(self.Parameters['enableFinal'])
            self.Drop_direction_Var.set(self.Parameters['direction'])
            self.paramDic['Scan Rate (V/s)_var'].set(self.Parameters['SR'])
            self.paramDic['Segments_var'].set(self.Parameters['seg'])
            self.Drop_eType_Var.set(self.Parameters['eType'])
            self.paramDic['Electrode Size (m)_var'].set(self.Parameters['eSize'])

            # Multi Load
            if self.MultiPlot_Active:
                self.Drop_Numplot_Var.set(self.MultiPlot['Numplot'])
                self.Drop_varType_Var.set(self.MultiPlot['varType'])
                self.Drop_varType2_Var.set(self.MultiPlot['varType2'])
                self.multiDic['Lowest Value_var'].set(self.MultiPlot['var_Low'])
                self.multiDic['Highest Value_var'].set(self.MultiPlot['var_High'])
                self.MultiPlot_Active = self.MultiPlot['Active']

            print('Load Done')

        except Exception as e:
            print('error loading save:', str(e))


########################################################################################################################
root = Tk()
root.option_add('*Font', FONT)
my_mclass = mclass(root)
root.mainloop()
