#This file is part of ChiVO, the Chilean Virtual Observatory
#A project sponsored by FONDEF (D11I1060)
#Copyright (C) 2015 Universidad Tecnica Federico Santa Maria Mauricio Solar
#                                                            Marcelo Mendoza
#                   Universidad de Chile                     Diego Mardones
#                   Pontificia Universidad Catolica          Karim Pichara
#                   Universidad de Concepcion                Ricardo Contreras
#                   Universidad de Santiago                  Victor Parada
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from Tkinter import *
import tkFileDialog
from ttk import *

import main as mn



class GUI(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   

        inputDir = ""
        outputDir = ""
        self.nuevo = Frame(self, relief=RAISED,borderwidth=1)
        self.E1 = Entry(self.nuevo, width=5)
        self.E2 = Entry(self.nuevo, width=5)
        self.E3 = Entry(self.nuevo, width=5)
        self.E4 = Entry(self.nuevo, width=5)
        self.visible = False
        self.parent = parent
        self.initUI()
        
    def initUI(self):
                  
        
        self.parent.title("Stacking V0.1b")
        self.style = Style()
        self.style.theme_use("default")
        
        iFrame = Frame(self, relief=RAISED, borderwidth=1)
        iFrame.pack(fill=BOTH, expand=1)
        oFrame = Frame(self, relief=RAISED, borderwidth=1)
        oFrame.pack(fill=BOTH, expand=1)

        self.pack(fill=BOTH, expand=1)      


        inputLabel = Label(iFrame, text="Input Directory:")
        inputLabel.pack(side=LEFT, padx=10)
        inputDirLabel = Label(iFrame, text="No selected")
        inputDirLabel.pack(side=LEFT, padx=5)
        inputButton = Button(iFrame, text="Choose", command=lambda:self.chooseDir(inputDirLabel,0))
        inputButton.pack(side=RIGHT, padx=5)

        outputLabel = Label(oFrame, text="Output Directory:")
        outputLabel.pack(side=LEFT, padx=5)
        outputDirLabel = Label(oFrame, text="No selected")
        outputDirLabel.pack(side=LEFT, padx=5)
        outputButton = Button(oFrame, text="Choose", command=lambda:self.chooseDir(outputDirLabel,1))
        outputButton.pack(side=RIGHT, padx=5)
        warningLabel = Label(self, text="Warning: The output directory's content will be deleted!")
        warningLabel.pack(side=LEFT, padx=5)

        closeButton = Button(self, text="Exit", command=self.close)
        closeButton.pack(side=RIGHT, padx=5)
        okButton = Button(self, text="Stack!", command=lambda:mn.stack(self.inputDir, self.outputDir,self.E1.get(),\
            self.E2.get(),self.E3.get(),self.E4.get()))
        okButton.pack(side=RIGHT)

        manualValues = Button(self, text = "Manual")
        manualValues.bind('<Button-1>', self.hideShow)
        manualValues.pack(side=RIGHT)
        


    def checkUncheck(self):
        print "variable is {0}".format(self.checked.get())

    def hideShow(self,event):
        if self.visible:
            self.nuevo.destroy()
            self.visible = False
        else:
            self.visible = True 
            self.nuevo = Frame(self.parent, relief=RAISED,borderwidth=1)
            self.nuevo.pack(fill=BOTH, expand=1)
            x1 = Label(self.nuevo, text="x1")
            x1.pack( side = LEFT,padx=20)
            self.E1 = Entry(self.nuevo, width=5)
            self.E1.pack(side = LEFT,padx=20)
            x2 = Label(self.nuevo, text="x2")
            x2.pack( side = LEFT,padx=20)
            self.E2 = Entry(self.nuevo, width=5)
            self.E2.pack(side = LEFT,padx=20)
            y1 = Label(self.nuevo, text="y1")
            y1.pack( side = LEFT,padx=20)
            self.E3 = Entry(self.nuevo, width=5)
            self.E3.pack(side = LEFT,padx=20)
            y2 = Label(self.nuevo, text="y2")
            y2.pack( side = LEFT,padx=20)
            self.E4 = Entry(self.nuevo, width=5)
            self.E4.pack(side = LEFT,padx=20)    

    def chooseDir(self,entry,inout):
        dirname = tkFileDialog.askdirectory(parent=self,initialdir=".",title='Please select a directory')
        if len(dirname ) > 0:
            print "You chose %s" % dirname 
            entry["text"] = dirname
            if inout == 0:
                self.inputDir = dirname
            else:
                self.outputDir = dirname
            return

    def close(self):
        self.parent.destroy()

def main():
    root = Tk()
    root.geometry("600x200+300+300")
    app = GUI(root)
    root.mainloop()  


if __name__ == '__main__':
    main()  
