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


import wx


def create(parent):
    return PySizes(parent)

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

    def __init__(self, parent):
        '''create the GUI and start the program'''
        self._init_ctrls(parent)


if __name__ == '__main__':
    app = wx.PySimpleApp()
    frame = create(None)
    frame.Show()

    app.MainLoop()
