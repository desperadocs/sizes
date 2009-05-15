#!/usr/bin/env python
#Boa:Frame:PySizes

'''GUI to run the sizes program for small-angle scattering analysis'''

########### SVN repository information ###################
# $Date$
# $Author$
# $Revision$
# $URL$
# $Id$
########### SVN repository information ###################


import sys, os, pprint, wx


def create(parent, cmdFile):
    return PySizes(parent, cmdFile)

#################################
###        Boa methods        ###
#################################

[wxID_PYSIZES, wxID_PYSIZESPANEL1, 
] = [wx.NewId() for _init_ctrls in range(2)]


class PySizes(wx.Frame):
    '''GUI to run the sizes program for small-angle scattering analysis'''

    def _init_ctrls(self, prnt):
        # generated method, don't edit
        wx.Frame.__init__(self, id=wxID_PYSIZES, name='PySizes', parent=prnt,
              pos=wx.Point(110, 145), size=wx.Size(400, 250),
              style=wx.DEFAULT_FRAME_STYLE, title='PySizes')
        self.SetClientSize(wx.Size(392, 216))

        self.panel1 = wx.Panel(id=wxID_PYSIZESPANEL1, name='panel1',
              parent=self, pos=wx.Point(0, 0), size=wx.Size(392, 216),
              style=wx.TAB_TRAVERSAL)

    def __init__(self, parent, cmdFile = None):
        '''create the GUI and start the program'''
        self._init_ctrls(parent)
        self.plotnum = 0
        self.params = self.__set_default_parameters__()
        self.__read_input_file__(cmdFile)
        pprint.pprint(self.params)
        #self.__write_input_file__('output.cmd')

    # ################################
    # ##       added methods       ###
    # ################################

    def __set_default_parameters__(self):
        '''set default parameters used by the program'''
        defaults = (
        "qMin       1.0e-8	      ",
        "qMax       100 	      ", #  1/A
        "contrast   1		      ",
        "dataFactor 1		      ", #  10^20/cm^4 (or 10^28/m^4)
        "esdFactor  1		      ", #    I = fac * ( i - bkg )
        "bkg        0.1 	      ", #  ESD = fac * err * esd
        "shape      1		      ", #    I = fac * ( i - bkg )
        "aspect     1.0 	      ", #  1=spheroids, no others yet
        "binType    1		      ", #  r x r x r*aspect
        "nRadii     100 	      ", #  1=Lin  0=Log
        "dMin       25		      ",
        "dMax       900 	      ", #  Angstroms
        "distType   1		      ",
        "defLevel   1.0e-8	      ", #  0=N(D)  1=f(D)  2=i(D)
        "iterMax    30		      ",
        "slitLength 0		      ", #  used only by MaxEnt
        "dE_E       0.0002	      ", #  for slit-smeared data
        "method     1		      ", #  incident wavelength spread
        "wdName     ./		      ", #  1=MaxEnt  0=regularization, 2=reg. with NNLS
        "project    test	      ",
        "sasFile    test.sas	      ",
        "cmdFile    test.cmd	      ",
        "font       *Times-Bold-R*14* ", # for plot titles
        )
        db = {}
        for pair in defaults:
            name, value = pair.split()
            db[name] = value
        return db

    def __read_input_file__(self, fileName):
        '''reads the parameters from an input command file'''
        fb = open(fileName, 'r')
        self.params['project'] = fb.readline().split()[0]
        self.params['sasFile'] = fb.readline().split()[0]
        self.params['qMin'], self.params['qMax'] = fb.readline().split()[0:2]
        self.params['contrast'] = fb.readline().split()[0]
        self.params['dataFactor'] = fb.readline().split()[0]
        self.params['esdFactor'] = fb.readline().split()[0]
        self.params['bkg'] = fb.readline().split()[0]
        self.params['shape'] = fb.readline().split()[0]
        self.params['aspect'] = fb.readline().split()[0]
        self.params['binType'] = fb.readline().split()[0]
        self.params['nRadii'] = fb.readline().split()[0]
        self.params['dMin'], self.params['dMax'] = fb.readline().split()[0:2]
        self.params['distType'] = fb.readline().split()[0]
        self.params['defLevel'] = fb.readline().split()[0]
        self.params['iterMax'] = fb.readline().split()[0]
        self.params['slitLength'] = fb.readline().split()[0]
        self.params['dE_E'] = fb.readline().split()[0]
        self.params['method'] = fb.readline().split()[0]
        fb.close()
        # example command file contents
        '''
             test-fm : Project Name  (only 1st item is read)
            test.sas : SAS file, contains columns: Q  i  esd
     1e-08       100 : qMin qMax, 1/A  (1.0e-8 to 100 means all data)
                 100 : rhosq       : scattering contrast, 10^20 1/cm^-4
                   1 : fac         :   I = fac * ( i - bkg )
                   3 : err         : ESD = fac * err * esd
                 0.1 : bkg         :   I = fac * ( i - bkg )
                   1 : shapeModel  (1=spheroids, no others yet)
                   1 : Aspect Ratio
                   0 : Bin Type    (1=Lin, 0=Log)
                  40 : nRadii
        25      9000 : dMin dMax, A
                   1 : n, in N(D)*V^n, 0=N(D), 1=f(D), 2=i(D)
              1.0e-6 : defaultDistLevel  (MaxEnt only)
                  91 : IterMax
                   0 : slitLength, 1/A
              0.0002 : dLambda/Lambda
                   1 : method (1=MaxEnt, 0=regularization)
        '''

    def __write_input_file__(self, fileName):
        '''writes the parameters to an input command file'''
        fmt1 = '%20s : %s\n'
        fmt2 = '%10s%10s : %s\n'
        fb = open(fileName, 'w')
        fb.write(fmt1 % (self.params['project'], 
            "Project Name, only 1st item is read"))
        fb.write(fmt1 % (self.params['sasFile'], 
            "SAS file, contains columns: Q  i  esd"))
        fb.write(fmt2 % (self.params['qMin'], self.params['qMax'], 
            "qMin qMax, 1/A  (1.0e-8 to 100 means all data)"))
        fb.write(fmt1 % (self.params['contrast'], 
            "rhosq       : scattering contrast, 10^20 1/cm^-4"))
        fb.write(fmt1 % (self.params['dataFactor'], 
            "fac         :   I = fac * ( i - bkg )"))
        fb.write(fmt1 % (self.params['esdFactor'], 
            "err         : ESD = fac * err * esd"))
        fb.write(fmt1 % (self.params['bkg'], 
            "bkg         :   I = fac * ( i - bkg )"))
        fb.write(fmt1 % (self.params['shape'], 
            "shapeModel  (1=spheroids, no others yet)"))
        fb.write(fmt1 % (self.params['aspect'], 
            "Aspect Ratio"))
        fb.write(fmt1 % (self.params['binType'], 
            "Bin Type    (1=Lin, 0=Log)"))
        fb.write(fmt1 % (self.params['nRadii'], 
            "nRadii"))
        fb.write(fmt2 % (self.params['dMin'], self.params['dMax'], 
            "dMin dMax, A"))
        fb.write(fmt1 % (self.params['distType'], 
            "n, in N(D)*V^n, 0=N(D), 1=f(D), 2=i(D)"))
        fb.write(fmt1 % (self.params['defLevel'], 
            "defaultDistLevel  (MaxEnt only)"))
        fb.write(fmt1 % (self.params['iterMax'], 
            "IterMax"))
        fb.write(fmt1 % (self.params['slitLength'], 
            "slitLength, 1/A"))
        fb.write(fmt1 % (self.params['dE_E'], 
            "dLambda/Lambda"))
        fb.write(fmt1 % (self.params['method'], 
            "method (0=reg., 1=MaxEnt, 2=reg+NNLS, 3=NNLS, 4=SVD)"))
        fb.close()

    # ################################
    # ##  event handling routines  ###
    # ################################



# ################################
# ##      other routines       ###
# ################################

def __usage__():
    '''suggest proper usage of this tool'''
    print "usage: pySize [command_file]"
    exit()


if __name__ == '__main__':
    if len(sys.argv) != 2:
        __usage__()
    cmdFile = sys.argv[1]
    app = wx.PySimpleApp()
    frame = create(None, cmdFile)
    frame.Show()

    app.MainLoop()
