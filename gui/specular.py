import wx
import time
import sys
import numpy
import matplotlib
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import wx.lib.inspection
#from wx.lib.mixins.listctrl import ColumnSorterMixin
from wx.lib.mixins.listctrl import ListCtrlAutoWidthMixin
import core.workspace as ws
from astropy.nddata import NDData
from astropy.table import Table
from matplotlib import pyplot as plt  

def load_test():
   ws.import_file("../fits/calibrated.ms.image.spectrum.J113740.6-010454.spw0.image.fits")
   #ws.import_file("../fits/calibrated.ms.line.spw0.source15.image.fits")
   #ws.import_file("../fits/NGC6240_continuum.fits")
#   ws.import_file("../fits/logfile_alma_hatlas_cycle1_inc-z_beye.fits")

class AutoWidthListCtrl(wx.ListCtrl, ListCtrlAutoWidthMixin):
    def __init__(self, parent):
        wx.ListCtrl.__init__(self, parent, -1, style=wx.LC_REPORT|wx.LC_ALIGN_LEFT|wx.LC_HRULES|wx.LC_VRULES)
        ListCtrlAutoWidthMixin.__init__(self)

class Explorer(wx.Frame):
   def __init__(self,parent,title):
      wx.Frame.__init__(self,parent,title=title,size=(800,500))
      self.createPanels()
      self.createMenuBar()
      self.createToolBar()

   def createToolBar(self):
      toolbar = self.CreateToolBar()
      qtool = toolbar.AddLabelTool(wx.ID_ANY, 'Quit', wx.Bitmap('quit.png'))
      toolbar.Realize()
      self.Bind(wx.EVT_TOOL, self.OnQuit, qtool)

   def createMenuBar(self):
      menubar = wx.MenuBar()
      fileMenu = wx.Menu()
      fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
      menubar.Append(fileMenu, '&File')
      self.SetMenuBar(menubar)
      self.Bind(wx.EVT_MENU, self.OnQuit, fitem)

   def updateWorkspace(self):
      items = ws.elements().items()
      self.myRowDict=dict()
      for key, data in items:
         #print key
         if isinstance(data,Table):
            index = self.wslist.InsertStringItem(sys.maxint, key)
            self.wslist.SetStringItem(index, 1, "Table")
            self.wslist.SetStringItem(index, 2, "---")
            self.myRowDict[index] = key
         elif  isinstance(data,NDData):
            (dim,shape,otype)=ws.real_dims(data)
            index = self.wslist.InsertStringItem(sys.maxint, key)
            self.wslist.SetStringItem(index, 1, otype)
            self.wslist.SetStringItem(index, 2, str(shape))
            self.myRowDict[index] = key

   def createPanels(self):
      
      self.figure = matplotlib.figure.Figure()
      self.mainPlot = self.figure.add_subplot(111)

      #Create main panels
      #self.mainPanel=wx.Panel(self)
      self.splitter = wx.SplitterWindow(self)
      self.innerPanel=wx.Panel(self.splitter)
      self.mainImgPanel=FigureCanvas(self.innerPanel,-1,self.figure)
      self.mainImgPanel.SetMinSize((10,10))
      x = np.arange(0, 5, 0.1);
      y1 = x
      y2 = 5 - x
      self.mainPlot.plot(y1,c="r")
      self.mainPlot.plot(y2,c="r")


      self.wslist = AutoWidthListCtrl(self.innerPanel)
      self.wslist.InsertColumn(0, 'Resource', width=130)
      self.wslist.InsertColumn(1, 'Type', width=70)
      self.wslist.InsertColumn(2, 'Dims', wx.LIST_FORMAT_RIGHT, 100)
      self.updateWorkspace()
      self.info = wx.TextCtrl(self.splitter, style=wx.TE_MULTILINE|wx.TE_READONLY,size=(500,400))    
      #Dynamics of the panel
      self.wslist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)

      vbox = wx.BoxSizer(wx.VERTICAL)
      vbox.Add(self.wslist,flag=wx.EXPAND|wx.TOP,border=5)
      vbox.Add(self.mainImgPanel,flag=wx.EXPAND,border=5,proportion=1)

      self.innerPanel.SetSizer(vbox)

      #Panel Layout
      #rbox = wx.BoxSizer(wx.HORIZONTAL)
      #rbox.Add(self.innerPanel,flag=wx.EXPAND,border=10)
      #rbox.Add(self.info,flag=wx.EXPAND,border=10)
      
      #self.mainPanel.SetSizer(rbox)
      
      self.splitter.SplitVertically(self.innerPanel, self.info)
      self.splitter.SetMinimumPaneSize(300)
      
      self.mainPlot.axis('off')
      self.mainPlot.axis('equal') 
      #Status Bar
      self.statusbar = self.CreateStatusBar()
      self.statusbar.SetStatusText("Ready.")

   def PlotMain(self):
      if isinstance(self.target,Table):
         x = np.arange(0, 5, 0.1);
         y1 = x
         y2 = 5 - x
         self.mainPlot.plot(y1,c="r")
         self.mainPlot.plot(y2,c="r")
         self.mainPlot.axis('off')
         self.mainPlot.axis('equal') 
         return
      ndd=self.target
      meta=self.target.meta
      (self.dim,self.shape,otype)=ws.real_dims(ndd)
      if self.dim == 1:
         self.mainPlot.plot(ndd.data)
         self.mainPlot.set_xlabel('index')
         self.mainPlot.set_ylabel('flux')
         self.mainPlot.axis('equal') 
      elif self.dim == 2:
         img=ndd.data
         if ndd.data.ndim==4:
            img=ndd.data[0][0]
         if ndd.data.ndim==3:
            img=img=ndd.data[0]
        
         dec_l=meta['CRVAL1'] - meta['CRPIX1']*meta['CDELT1']
         dec_u=meta['CRVAL1'] + (meta['NAXIS1']- meta['CRPIX1'])*meta['CDELT1']
         ra_l=meta['CRVAL2'] - meta['CRPIX2']*meta['CDELT2']
         ra_u=meta['CRVAL2'] + (meta['NAXIS2']- meta['CRPIX2'])*meta['CDELT2']
         self.mainPlot.imshow(img,extent=[dec_l,dec_u,ra_l,ra_u])
         self.mainPlot.set_xlabel(meta['CTYPE1'])
         self.mainPlot.set_ylabel(meta['CTYPE2'])
         #self.set_xticks(np.arange(max_val))
         #ax.set_yticks(np.arange(max_val))

      elif self.dim == 3:
         if ndd.data.ndim==4:
            img=ndd.data[0]
         else:
            img=ndd.data
         #self.mainPlot.imshow(img.mean(axis=(0)))
         dec_l=meta['CRVAL1'] - meta['CRPIX1']*meta['CDELT1']
         dec_u=meta['CRVAL1'] + (meta['NAXIS1']- meta['CRPIX1'])*meta['CDELT1']
         ra_l=meta['CRVAL2'] - meta['CRPIX2']*meta['CDELT2']
         ra_u=meta['CRVAL2'] + (meta['NAXIS2']- meta['CRPIX2'])*meta['CDELT2']
         self.mainPlot.imshow(img.mean(axis=(0)),extent=[dec_l,dec_u,ra_l,ra_u])
         self.mainPlot.set_xlabel(self.target.meta['CTYPE1'])
         self.mainPlot.set_ylabel(self.target.meta['CTYPE2'])
      self.figure.tight_layout()
 
   def getNiceDictRepr(self,aDict):
      return '\n'.join('%s = %s' % t for t in aDict.iteritems())

   def onItemSelected(self,event):
      currentItem = event.m_itemIndex
      key=self.myRowDict[currentItem]
      self.statusbar.SetStatusText("Loading Preview...")
      items=ws.elements()
      self.mainPlot.clear()
      self.target=items[key]
      self.PlotMain()
      self.mainImgPanel.draw()
      self.mainImgPanel.Refresh()
      self.wslist.Refresh()
      self.info.Clear()
      self.info.WriteText(self.getNiceDictRepr(self.target.meta))
      self.info.SetInsertionPoint(0)
      self.info.Refresh()
      self.statusbar.SetStatusText("Ready.")

   def OnQuit(self, e):
      self.Close()

   def Load(self,event):
     self.statusbar.SetStatusText("Loading...")
     time.sleep(3)
     self.statusbar.SetStatusText("Ready.")

#class Specular(wx.Frame):
#   def __init__(self,parent,title):
#      wx.Frame.__init__(self,parent,title=title,size=(1024,768))
#      self.createPanels()
#      self.createMenuBar()
#      self.createToolBar()
#
#   def createToolBar(self):
#      toolbar = self.CreateToolBar()
#      qtool = toolbar.AddLabelTool(wx.ID_ANY, 'Quit', wx.Bitmap('quit.png'))
#      toolbar.Realize()
#      self.Bind(wx.EVT_TOOL, self.OnQuit, qtool)
#
#   def createMenuBar(self):
#      menubar = wx.MenuBar()
#      fileMenu = wx.Menu()
#      fitem = fileMenu.Append(wx.ID_EXIT, 'Quit', 'Quit application')
#      menubar.Append(fileMenu, '&File')
#      self.SetMenuBar(menubar)
#      self.Bind(wx.EVT_MENU, self.OnQuit, fitem)
#
#   def createFigures(self):
#      self.figure1 = matplotlib.figure.Figure()
#      self.figure2 = matplotlib.figure.Figure()
#      self.figure3 = matplotlib.figure.Figure()
#      self.mainPlot = self.figure1.add_subplot(111)
#      self.auxPlot  = self.figure2.add_subplot(111)
#      self.resultPlot  = self.figure3.add_subplot(111)
#   
#   def selectPixel(self, event):
#        print 'in selectPixel callback: clicked at (%g, %g)' % (event.x, event.y)
#        #menu = wxMenu()
#        #item_id = wxNewId()
#        #menu.Append(item_id, 'item') 
#        #wx.EVT_MENU(menu, item_id, self.callback)
#        #self.PopupMenu(menu, wx.Point(event.x, event.y))
#        #menu.Destroy()
#
#   def createPanels(self):
#      self.createFigures()
#      
#      #Create main panels
#      self.primary_panel = wx.Panel(self)
#      main_box = wx.BoxSizer(wx.HORIZONTAL)
#      self.mainPanel=wx.Panel(self)
#      inner_box = wx.BoxSizer(wx.VERTICAL)
#     
#      stc=wx.StaticBox(self.mainPanel,label="Main View")
#      self.mainImgPanel=FigureCanvas(self.mainPanel,-1,self.figure1)
#      self.mainImgPanel.SetMinSize((10,10))
#      self.mainImgPanel.mpl_connect('button_press_event', self.selectPixel) 
#      self.sld = wx.Slider(self.mainPanel, value=50, minValue=0, maxValue=100,style=wx.SL_HORIZONTAL|wx.SL_LABELS)
#      self.sld.Enable(False)
#      self.sld.Bind(wx.EVT_SCROLL, self.OnSliderScroll)
#      inner_box.Add(stc)
#      inner_box.Add(self.mainImgPanel,flag=wx.EXPAND|wx.TOP,proportion=1)
#      inner_box.Add(self.sld,flag=wx.EXPAND|wx.BOTTOM)
#      self.mainPanel.SetSizer(inner_box)
#
#      self.auxPanel=FigureCanvas(self.primary_panel,-1,self.figure2)
#      self.auxPanel.SetMinSize((10,10))
#
#      #Create list of objects
#      self.workspace = wx.Panel(self)
#      inner_box = wx.BoxSizer(wx.VERTICAL)
#      
#      self.wslist = AutoWidthListCtrl(self.primary_panel)
#      self.wslist.InsertColumn(0, 'Resource')#, width=500)
#      self.wslist.InsertColumn(1, 'Type')#, width=40)
#      self.wslist.InsertColumn(2, 'Dims')#, wx.LIST_FORMAT_RIGHT, 30)
#      self.updateWorkspace()
#      
#      #Dynamics of the panel
#      self.wslist.Bind(wx.EVT_LIST_ITEM_SELECTED, self.onItemSelected)
#
#      #Panel Layout
#      main_box.Add(self.wslist,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
#      main_box.Add(self.mainPanel,flag=wx.EXPAND,proportion=1,border=10)
#      main_box.Add(self.auxPanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
#      self.primary_panel.SetSizer(main_box)
#      
#      #Create analysis panels
#      self.secondary_panel = wx.Panel(self, 1)
#      second_box = wx.BoxSizer(wx.HORIZONTAL)
#      self.analysisPanel=AutoWidthListCtrl(self.primary_panel)
#      self.resultPanel=FigureCanvas(self.secondary_panel,-1,self.figure3)
#      self.resultPanel.SetMinSize((10,10))
#      
#      #Panel Layout
#      second_box.Add(self.info,flag=wx.EXPAND|wx.LEFT,proportion=1,border=10)
#      second_box.Add(self.analysisPanel,flag=wx.EXPAND,proportion=1,border=10)
#      second_box.Add(self.resultPanel,flag=wx.EXPAND|wx.RIGHT,proportion=1,border=10)
#      self.secondary_panel.SetSizer(second_box)
#
#      #General Layout
#      frame_box = wx.BoxSizer(wx.VERTICAL)
#      frame_box.Add(self.primary_panel,flag= wx.EXPAND|wx.TOP,proportion=1,border=10)
#      frame_box.Add(self.secondary_panel,flag= wx.EXPAND|wx.BOTTOM,proportion=1,border=10)
#      self.SetSizer(frame_box)
#      
#      #Status Bar
#      self.statusbar = self.CreateStatusBar()
#      self.statusbar.SetStatusText("Ready.")
#
#   def PlotMain(self,ndd):
#      (self.dim,self.shape,otype)=ws.real_dims(ndd)
#      self.sld.SetValue(50)
#      self.sld.SetMax(100)
#      if self.dim == 1:
#         self.mainPlot.plot(ndd)
#         self.sld.Enable(False)
#      elif self.dim == 2:
#         img=ndd
#         if ndd.ndim==4:
#            img=ndd[0][0]
#         if ndd.ndim==3:
#            img=ndd[0]
#         self.mainPlot.imshow(img)
#         
#         self.sld.Enable(False)
#      elif self.dim == 3:
#         if ndd.ndim==4:
#            self.target=ndd[0]
#         self.sld.SetMax(self.shape[0]-1)
#         self.sldIndex=int(self.shape[0]/2)
#         self.sld.SetValue(self.sldIndex)
#         self.sld.Enable(True)
#         self.mainPlot.imshow(self.target[self.sldIndex])
#         self.meanFlux=self.target.data.mean(axis=(1,2))
#         self.auxPlot.plot(self.meanFlux)
#         x1,x2,y1,y2 = self.auxPlot.axis()
#         self.auxPlot.axvline(self.sldIndex,c="m")
#         self.auxPlot.axis((x1,x2,self.target.data.min(),self.target.data.max()))
# 
#   def OnSliderScroll(self, e):
#      obj = e.GetEventObject()
#      self.sldIndex = obj.GetValue()
#      self.mainPlot.clear()
#      self.auxPlot.clear()
#      self.mainPlot.imshow(self.target[self.sldIndex])
#      self.auxPlot.plot(self.meanFlux)
#      x1,x2,y1,y2 = self.auxPlot.axis()
#      self.auxPlot.axis((x1,x2,self.target.data.min(),self.target.data.max()))
#      self.auxPlot.axvline(self.sldIndex,c="r")
#      self.mainImgPanel.draw()
#      self.mainImgPanel.Refresh()
#      self.auxPanel.draw()
#      self.auxPanel.Refresh()


load_test()
app = wx.App(redirect=False)
frame=Explorer(None,"Specular Explorer")
frame.Show()
frame.Centre()
wx.lib.inspection.InspectionTool().Show()
app.MainLoop()
