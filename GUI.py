import matplotlib as plt
plt.use("TkAgg") 
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure

import tkinter as tk
from tkinter import Tk, ttk

LARGE_FONT=("Verdana",12)

class Gui(Tk):
    """
    Initializing the GUI
    """
    def __init__(self, *args, **kwargs) -> None:
        Tk.__init__(self, *args, **kwargs)

        Tk.iconbitmap(self,default='')  #get the icon and move into folder
        Tk.wm_title(self,"Drillsoft Multiphase Flow Temperature Simulation")

        container=tk.Frame(self)
        container.pack(side='top',fill='both',expand=True)
        container.grid_rowconfigure(0,weight=1)
        container.grid_columnconfigure(0,weight=1)

        self.frames={}

        for F in (StartPage,PageOne,PageTwo,PageThree):
            frame=F(container, self)
            self.frames[F]=frame
            frame.grid(row=0,column=0,sticky='nsew')

        self.show_frame(StartPage)

    def show_frame(self,cont): #cont=controller
        frame=self.frames[cont]
        frame.tkraise()

class StartPage(tk.Frame):
    """
    Pages are created and stored here
    """
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        label=tk.Label(self,text='Start Page',font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1=ttk.Button(self,text="Visit Page 1",
                            command=lambda: controller.show_frame(PageOne)) #lambda here is important
        button1.pack()

        button2=ttk.Button(self,text="Page Two",
                            command=lambda: controller.show_frame(PageTwo))
        button2.pack()

        button3=ttk.Button(self,text="Page Three",
                            command=lambda: controller.show_frame(PageThree))
        button3.pack()


class PageOne(tk.Frame):

    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        label=tk.Label(self,text="Page One",font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1=ttk.Button(self,text="back to home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2=ttk.Button(self,text="Page Two",
                            command=lambda: controller.show_frame(PageTwo))
        button2.pack()

class PageTwo(tk.Frame):

    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        label=tk.Label(self,text="Page Two",font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1=ttk.Button(self,text="back to home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        button2=ttk.Button(self,text="Page One",
                            command=lambda: controller.show_frame(PageOne))
        button2.pack()

class PageThree(tk.Frame):

    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        label=tk.Label(self,text="Graph Page",font=LARGE_FONT)
        label.pack(pady=10,padx=10)

        button1=ttk.Button(self,text="back to home",
                            command=lambda: controller.show_frame(StartPage))
        button1.pack()

        f=Figure(figsize=(5,5),dpi=100)
        a=f.add_subplot(111)
        a.plot([1,2,3,4,5],[4,5,1,2,3])

        canvas=FigureCanvasTkAgg(f,self)
        canvas.draw()

        toolbar=NavigationToolbar2Tk(canvas,self)
        toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP,fill=tk.BOTH,expand=True)
        canvas.get_tk_widget().pack(side=tk.TOP,fill=tk.BOTH,expand=True)


app=Gui()
app.mainloop()

