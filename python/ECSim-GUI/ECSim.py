import pyecsim as ecs
import matplotlib.pyplot as plt
import numpy as np
import math as m

DataDic = {}

def ECSim(MultiPlot, Parameters, Electron_Step, Chemical_Step, Count):


    sim = ecs.Simulation(True) # or False to suppress output

    # Electron Steps
    for rdx in Electron_Step.keys():
        sim.sys.addRedox(Electron_Step[rdx])

    # Chemical Steps
    for rxn in Chemical_Step.keys():
        sim.sys.addReaction(Chemical_Step[rxn])

    # Electrode Settings
    sim.el.disk(3.0e-3) # radius [m]

    # CV settings
    sim.exper.setScanPotentials(Parameters['init_E'], Parameters['segList'], Parameters['final_E']) # potentials [V]: initial, [0 or more vertices], final
    sim.exper.setScanRate(Parameters['SR']) # scan rate [V/s]

    # Plot Settings
    [potential, current] = sim.run()
    MultiPlot['Plots']['plt'+str(len(MultiPlot['Plots'])+1)] = [potential, current, MultiPlot['Label']]
    plt.plot(potential, current, label=MultiPlot['Label'])
    DataDic['potential' + str(Count)] = potential
    DataDic['current' + str(Count)] = current
    return [MultiPlot, DataDic]


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

def RunSim(MultiPlot, Parameters, Electron_Step, Chemical_Step):

###################
# Segment Creator #
###################

    for i in range(Parameters['seg']):
        if Parameters['direction'] == 'pos':
            if(i % 2 == 0):
                Parameters['segList'].append(Parameters['high_E'])
            else:
                Parameters['segList'].append(Parameters['low_E'])

        elif Parameters['direction'] == 'neg':
            if (i % 2 == 0):
                Parameters['segList'].append(Parameters['low_E'])
            else:
                Parameters['segList'].append(Parameters['high_E'])
    if Parameters['enableFinal'] == True:
        del Parameters['segList'][-1]

    if Parameters['enableFinal'] == False:
        Parameters['final_E'] = Parameters['segList'][-1]

    print(Parameters['segList'])


#########################################
# Log/Normal Distribution of Parameters #
#########################################

    if MultiPlot['Active']:

        if any(varType in MultiPlot['varType'] for varType in ['k_f', 'k_b', 'alpha', 'k_e', 'SR']):
            MultiPlot['var_Low'] = m.log(MultiPlot['var_Low'], 10)
            MultiPlot['var_High'] = m.log(MultiPlot['var_High'], 10)

        varRange = np.linspace(MultiPlot['var_Low'], MultiPlot['var_High'], MultiPlot['Numplot'])

        print(varRange,'This is before 10^x')


        if any(varType in MultiPlot['varType'] for varType in ['k_f', 'k_b', 'alpha', 'k_e', 'SR']):
            varRange = [10**x for x in varRange]

        print(varRange, 'This is after 10^x')


#################################
# MultiPlot Variable Assignment #
#################################

        for n in range(MultiPlot['Numplot']):

            if MultiPlot['varType'] in Parameters.keys():
                Parameters[MultiPlot['varType']] = varRange[n]

            elif MultiPlot['varType'] == 'k_f':
                Chemical_Step[MultiPlot['varType2']].setKf(varRange[n])

            elif MultiPlot['varType'] == 'k_b':
                Chemical_Step[MultiPlot['varType2']].setKb(varRange[n])


            elif MultiPlot['varType'] == 'k_e':
                Electron_Step[MultiPlot['varType2']].setKe(varRange[n])

            elif MultiPlot['varType'] == 'alpha':
                Electron_Step[MultiPlot['varType2']].setAlpha(varRange[n])

        #############################
        # Plot and Label Assignment #
        #############################

            MultiPlot['Label'] = ''.join([
                ''.join(MultiPlot['varType']),'=', '{:.2e}'.format(varRange[n]),
                # 'k_e=', '{:.2e}'.format(vars(Electron_Step[MultiPlot['varType'][0]])['k_e']), ', ',
                # 'alpha=', '{:.2e}'.format(vars(Electron_Step[MultiPlot['varType'][0]])['alpha']), ', ',
            ])
            Data = ECSim(MultiPlot, Parameters, Electron_Step, Chemical_Step, n)
            MultiPlot = Data[0]
            DataDic.update(Data[1])
            print(MultiPlot['Plots'].keys())
        # plt.text(0.02,0.5,[Parameters, vars(Solution['ox1'])], fontsize=14, transform=plt.gcf().transFigure)
        # plt.xlabel("Potential (V)")
        # plt.ylabel("Current (A)")
        # plt.legend(loc="upper left")
        # plt.grid('on', linestyle='--')
        # plt.show()
        return [MultiPlot['Plots'], DataDic]
    else:
        MultiPlot['Label'] = ''.join([
            'Plot #1: ',
            'SR=', str(Parameters['SR']), ', ',
        ])
        Data = ECSim(MultiPlot, Parameters, Electron_Step, Chemical_Step, 0)
        MultiPlot = Data[0]
        DataDic.update(Data[1])
        print(MultiPlot['Plots'].keys())
        return [MultiPlot['Plots'], DataDic]

########################################################################################################################
########################################################################################################################
########################################################################################################################


def VarSend(Multi, Param, Elec, Chem):

    Plots = RunSim(Multi, Param, Elec, Chem)
    return Plots


"""
    MultiPlot = {
        'Numplot': 1,
        'varType': ['rdx1','alpha'],
        'var_Low': 0.5,
        'var_High': 0.6,
        'Label': '',
        'varValue': 0,
        'varTest': 0,
    }

    Parameters = {
        'init_E': -1,
        'high_E': 1.7,
        'low_E': -1.0,
        'final_E': -1,
        'enableFinal': False,
        'direction': 'pos',
        'SR': 0.00001,
        'seg': 2,
        'segList': [],
        'eType': 'disk',
        'eSize': 3.0e-3,
    }

    # name, concentration [mol/m3], diffusion coefficient [m2/s]
    Solution = {
        'ox2-OH2': ecs.Species('Ru(II)-OH2', 0.05, 7.29e-10),
        'ox3-OH2': ecs.Species('Ru(III)-OH2', 0.0, 7.29e-10),  # 1st e step
        'ox3-OH': ecs.Species('Ru(III)-OH', 0.0, 7.29e-10),  # Rxn ox3-OH2 -> ox3-OH + prod1 (H+)
        'ox4-OH': ecs.Species('Ru(IV)-OH', 0.0, 7.29e-10),  # 2nd e step
        'ox4-O': ecs.Species('Ru(IV)-O', 0.0, 7.29e-10),  # Rxn ox4-OH -> ox4-O + H+
        'ox5-O': ecs.Species('Ru(V)-O', 0.0, 7.29e-10),  # 3rd e step
        'ox3-O-OH': ecs.Species('Ru(III)-O-OH', 0.0, 7.29e-10),  # Rxn ox5-O -> ox3-O-OH + H+
        'ox4-O-OH': ecs.Species('Ru(IV)-O-OH', 0.0, 7.29e-10),  #Rxn ox4-O-OH -> ox4-O-O + H+
        'ox4-O-O': ecs.Species('Ru(IV)-O-O', 0.0, 7.29e-10),  # Rxn ox4-O-O -> O2 + ox2
        'ox2': ecs.Species('Ru(II)', 0.0, 7.29e-10),  # ox2 -> ox3, 4th e step
        'ox3': ecs.Species('Ru(IV)-O-O', 0.0, 7.29e-10),  # Rxn ox3 + H2O -> ox3-OH2
        'O2': ecs.Species('Ru(IV)-O-O', 0.0, 2.0e-9),  # Googled this
        # 'H+': ecs.Species('H+', 1e-7, 7.62e-9),  # D for H+ found here: http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
    }

    # ox, red, n, E_0 [V], k_e [m/s], alpha
    Electron_Step = {
        'rdx1': ecs.Redox(Solution['ox3-OH2'], Solution['ox2-OH2'], 1, 0.5, 0.01, 0.5).enable(),  # Ru(II)-OH2 -> Ru(III)-OH2
        'rdx2': ecs.Redox(Solution['ox4-OH'], Solution['ox3-OH'], 1, 1.5, 0.01, 0.5).enable(),  # Ru(III)-OH -> Ru(IV)-OH
        'rdx3': ecs.Redox(Solution['ox5-O'], Solution['ox4-O'], 1, 1.4, 0.01, 0.5).enable(),  # Ru(IV)-O -> Ru(V)-O
        'rdx4': ecs.Redox(Solution['ox4-O-OH'], Solution['ox3-O-OH'], 1, 1.2, 100000, 0.5).enable(),  # Ru(III)-O-OH -> Ru(IV)-O-OH
    }

    # reactant1, reactant2, product1, product2, k_f, k_b
    Chemical_Step = {
        'rxn1': ecs.Reaction(Solution['ox3-OH2'], None, Solution['ox3-OH'], None, 50, 5).enable(),
        'rxn2': ecs.Reaction(Solution['ox4-OH'], None, Solution['ox4-O'], None, 100000, 10).enable(),
        'rxn3': ecs.Reaction(Solution['ox5-O'], None, Solution['ox3-O-OH'], None, 100000, 10).enable(),
        'rxn4': ecs.Reaction(Solution['ox4-O-OH'], None, Solution['ox4-O-O'], None, 100000, 10).enable(),
        'rxn5': ecs.Reaction(Solution['ox4-O-O'], None, Solution['ox2'], Solution['O2'], 1000000000, 0).enable(),
        'rxn6': ecs.Reaction(Solution['ox2'], None, Solution['ox2-OH2'], None, 10000000, 0).enable(),
    }
"""
if __name__ == '__main__':

    MultiPlot = {
        'Numplot': 3,
        'varType': ['SR'],
        'var_Low': 0.01,
        'var_High': 0.0001,
        'Label': '',
        'varValue': 0,
        'varTest': 0,
    }

    Parameters = {
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

    # name, concentration [mol/m3], diffusion coefficient [m2/s]
    Solution = {
        'ox2-OH2': ecs.Species('Ru(II)-OH2', 0.05, 7.29e-10),
        'ox3-OH2': ecs.Species('Ru(III)-OH2', 0.0, 7.29e-10),  # 1st e step
        'ox3-OH': ecs.Species('Ru(III)-OH', 0.0, 7.29e-10),  # Rxn ox3-OH2 -> ox3-OH + prod1 (H+)
        'ox4-OH': ecs.Species('Ru(IV)-OH', 0.0, 7.29e-10),  # 2nd e step
        'ox4-O': ecs.Species('Ru(IV)-O', 0.0, 7.29e-10),  # Rxn ox4-OH -> ox4-O + H+
        'ox5-O': ecs.Species('Ru(V)-O', 0.0, 7.29e-10),  # 3rd e step
        'ox3-O-OH': ecs.Species('Ru(III)-O-OH', 0.0, 7.29e-10),  # Rxn ox5-O -> ox3-O-OH + H+
        'ox4-O-OH': ecs.Species('Ru(IV)-O-OH', 0.0, 7.29e-10),  #Rxn ox4-O-OH -> ox4-O-O + H+
        'ox4-O-O': ecs.Species('Ru(IV)-O-O', 0.0, 7.29e-10),  # Rxn ox4-O-O -> O2 + ox2
        'ox2': ecs.Species('Ru(II)', 0.0, 7.29e-10),  # ox2 -> ox3, 4th e step
        'ox3': ecs.Species('Ru(III)', 0.0, 7.29e-10),  # Rxn ox3 + H2O -> ox3-OH2
        'O2': ecs.Species('Ru(IV)-O-O', 0.0, 2.0e-9),  # Googled this
        # 'H+': ecs.Species('H+', 1e-7, 7.62e-9),  # D for H+ found here: http://omh.umeche.maine.edu/pdfs/JChemPhys_135_124505.01pdf.pdf
    }

    # ox, red, n, E_0 [V], k_e [m/s], alpha
    Electron_Step = {
        'rdx1': ecs.Redox(Solution['ox3-OH2'], Solution['ox2-OH2'], 1, 0.55, 0.1, 0.5).enable(),  # Ru(II)-OH2 -> Ru(III)-OH2
        'rdx2': ecs.Redox(Solution['ox4-OH'], Solution['ox3-OH'], 1, 1.3, 0.1, 0.5).enable(),  # Ru(III)-OH -> Ru(IV)-OH
        'rdx3': ecs.Redox(Solution['ox5-O'], Solution['ox4-O'], 1, 1.3, 0.1, 0.5).enable(),  # Ru(IV)-O -> Ru(V)-O
        'rdx4': ecs.Redox(Solution['ox4-O-OH'], Solution['ox3-O-OH'], 1, 1.2, 0.1, 0.5).enable(),  # Ru(III)-O-OH -> Ru(IV)-O-OH
        'rdx5': ecs.Redox(Solution['ox3'], Solution['ox2'], 1, 0.2, 0.1, 0.5).enable(),
    }
    s = 1e3
    r1 = 1e3      # Ru(III)-OH2 -> Ru(III)-OH
    r2 = 1e-1     # Ru(IV)-OH -> Ru(IV)=O
    r3 = 1e2
    r4 = 1e2
    r5 = 1e3
    r6 = 1e-1     # Ru(II) -> Ru(II)-OH2
    r7 = 1e8      # Ru(III) -> Ru(III)-OH2

    #     # reactant1, reactant2, product1, product2, k_f, k_b
    Chemical_Step = {
        'rxn1': ecs.Reaction(Solution['ox3-OH2'], None, Solution['ox3-OH'], None, 1*r1*s, 0.1*r1*s).enable(),
        'rxn2': ecs.Reaction(Solution['ox4-OH'], None, Solution['ox4-O'], None, 1*r2*s, 0.1*r2*s).enable(),
        'rxn3': ecs.Reaction(Solution['ox5-O'], None, Solution['ox3-O-OH'], None, 1*r3*s, 0*r3*s).enable(),
        'rxn4': ecs.Reaction(Solution['ox4-O-OH'], None, Solution['ox4-O-O'], None, 1*r4*s, 1*r4*s).enable(),
        'rxn5': ecs.Reaction(Solution['ox4-O-O'], None, Solution['ox2'], Solution['O2'], 1*r5*s, 0*r5*s).enable(),
        'rxn6': ecs.Reaction(Solution['ox2'], None, Solution['ox2-OH2'], None, 1*r6*s, 0.1*r6*s).enable(),  # Make slow kf > kb
        'rxn7': ecs.Reaction(Solution['ox3'], None, Solution['ox3-OH2'], None, 1*r7*s, 0.1*r7*s).enable(),  # Make fast kf > kb

    }

    RunSim(MultiPlot, Parameters, Solution, Electron_Step, Chemical_Step)
