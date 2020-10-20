import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
import turtle as t1
import time
import math
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import tkinter as tk

window = tk.Tk()
window.title('Nanojunction Builder')

window.geometry("675x900+0+0")
window.resizable(0, 0) #Don't allow resizing in the x or y direction
wn = t1.Screen()
wn.setup(600,900,startx=680, starty=0)
wn.title("Nanojunction Builder")

border_effects = {

    "flat": tk.FLAT,

    "sunken": tk.SUNKEN,

    "raised": tk.RAISED,

    "groove": tk.GROOVE,

    "ridge": tk.RIDGE,

}
texts =["Import Slab1","Import Slab2", "Import XYZ", "Export XYZ","Bond Slab1","Bond Slab2","Coincidence Lat.","Clear Coinc. Lat.", "Render Ovito","Edges", "Surface","Volume", "Interaction",
        "Select Pivots","Del Overlap","Create Bonds"]

coordx_L = []
coordy_L = []
coordz_L = []
coordx_R = []
coordy_R = []
coordz_R = []
coordx_L_trans = []
coordy_L_trans = []
coordz_L_trans = []
coordx_R_trans = []
coordy_R_trans = []
coordz_R_trans = []
coordpx_L = []
coordpy_L = []
coordpz_L = []

coordpx_R = []
coordpy_R = []
coordpz_R = []

atoms_sel = []
atomsT_L = []
atomsT_R=[]
atoms_L = []
atoms_R = []
atomsT =[]
new_Lattice = []
xwin1 = -250
xwin2 = 250
ywin1 = -250
ywin2 = 250
alpha = 10
x1 = 0.0
y1 = 0.0
x2= 0.0
y2 = 0.0
x1_y1_flag = False
beta_L=[0,0.2]
beta_R=[0,0.2]
xPix1 = 0.0 
yPix1 = 0.0
    
    
xPix2 = 0
yPix2 = 0
dist = 0.0
angle = 0
x_click = 0.0
y_click = 0.0
direction = " "
bondsT_L = []
bondsT_R = []
bondsT= []
midpointsx= []
midpointsy= []
num_slab = 2
slab_status = [False]*num_slab
bond_status = [False]*num_slab
reset = False

bondT = t1.Turtle()
    
bondT.ht()
bondT.penup()
bondT.pensize(1)
bondT.speed(0)
bondT.ht()
bondT1 = bondT.clone()
bondT1.ht()
bondT2 = bondT.clone()
bondT2.ht()
class Shape(t1.Turtle):
    def __init__(self,shape):
        
        t1.Turtle.__init__(self)
        self.ht()
        self.clear()
        self.shape(shape)
        
        self.penup()
        self.speed(0)
        self.shapesize(0.3,0.3)
        
    def del_atom(self):
        self.ht()
        atomsT.remove(self)       
    
    def dragging(self,x,y):  
        self.ondrag(None)
        self.setheading(self.towards(x, y))
        self.goto(x, y)
        self.ondrag(dragging)

    def clickRight(x,y):
        self.clear()

# used for atoms        
atomT = Shape("circle")

# used for coincidence lattice points
atomC = Shape("square")
atomC.ht()
atomC.shapesize(0.6,0.6)
atomC.color("black")

atomC.fillcolor("")

atomC.pensize(4)
            
# Convert an (x,y) coordinate from Angstrom to pixels
def convert_angstrom_pixels(vx,vy,alpha,beta): #,betax,betay
    #global xwin1,xwin2,ywin1,ywin2
    global x1,x2,y1,y2
    global xPix1,yPix1,xPix2,yPix2,x1_y1_flag
    global slab_status
    #sets the size of the atomic system in angstrom xmin,ymin... xmax,ymax
   
    coordx = vx
    coordy = vy
    
    pixelx = []
    pixely = []
    
    if x1_y1_flag == False:
        
        x1 = min(coordx)
        y1 = min(coordy)
        x2 = max(coordx)
        y2 = max(coordy)   
        x1_y1_flag =True   
    #print("x1=",x1,"x2=",x2)
    
    xPix1 = xwin1 + beta[0] *(xwin2-xwin1)  
    yPix1 = ywin1 + beta[1] *(ywin2-ywin1)
    
    
    xPix2 = xPix1 + alpha *(x2-x1)
    yPix2 = yPix1 + alpha *(y2-y1)
      
    for i in range(len(coordx)): 
        xPix= xPix1+((xPix2- xPix1)/(x2-x1))*(coordx[i]-x1)
        yPix= yPix1+((yPix2- yPix1)/(y2-y1))*(coordy[i]-y1)
           
        
        pixelx.append(xPix)
        pixely.append(yPix)
            
    return pixelx,pixely     


def OpenFile():
    filename = askopenfilename()
    print(filename)
    return filename

def SavetoFile():
    
    filename = asksaveasfilename(defaultextension='.txt')
    print(filename)
    return filename

#Export an atomic structure
def export_file():
    filename = SavetoFile()
    file =open(filename,"w")
    for i in range(len(coordx_L)):
        file.write(str(coordx_L[i])+","+ str(coordy_L[i])+","+ str(coordz_L[i])+","+str(i+1)+","+"L"+"\n")
    for i in range(len(coordx_R)):
        file.write(str(coordx_R[i])+","+ str(coordy_R[i])+","+ str(coordz_R[i])+","+str(i+1)+","+"R"+"\n")    
    file.close()

# Import a Slab   
def plot_sheets(alpha,beta,slab): #betax,betay,
    #import turtle as t
    #global frame1,canvas1,t1,wn
    global coordx_L,coordy_L,coordz_L
    global coordx_R,coordy_R,coordz_R
    global coordpx_L,coordpy_L,coordpz_L
    global coordpx_R,coordpy_R,coordpz_R
    global atomsT_L,atomsT_R
    global atoms_L,atoms_R, slab_status,atomT,reset

    try:
 
        #read a file x,y,z
        filename = OpenFile()
        file =open(filename,"r")
        
        
        if slab == 1:

            atomT.color("blue")
            slab_status[0] =True
        if slab ==2:    

            atomT.color("red")
            atomT.fillcolor("")
            atomT.shapesize(0.3,0.3)
            slab_status[1] =True
    except:
        print("file does not exist")
        return None
   
    for rec in file:
        
        line = rec.split(",")
        a=line[len(line)-1].splitlines()
        line[len(line)-1] =a[0]
        
        x =  float(line[0]) 
        y =  float(line[1])
        z =  float(line[2])
        #add atom number, and L,R may be all a dictionary
        if slab == 1:
           coordx_L.append(x)
           coordy_L.append(y)
           #print(x,y)
           coordz_L.append(z)
           
        if slab == 2:
           coordx_R.append(x)
           coordy_R.append(y)
           #print(x,y)
           coordz_R.append(z)    
    file.close()
    # conversion from Angstrom to pixels
    if slab == 1:
        coordpx_L,coordpy_L=convert_angstrom_pixels(coordx_L,coordy_L,alpha,beta)#betax,betay
        
    if slab == 2:
        coordpx_R,coordpy_R=convert_angstrom_pixels(coordx_R,coordy_R,alpha,beta)#betax,betay

        
    if slab == 1:
        coordpx  = coordpx_L
        coordpy  = coordpy_L
   
    if slab == 2:
        coordpx  = coordpx_R
        coordpy  = coordpy_R
    wn.tracer(0)


    for i in range(len(coordpx)):
        
        atomT1 = atomT.clone()
        
        #atomT1.ht()
        atomT1.penup()
        atomT1.speed(0)
        atomT1.goto(coordpx[i],coordpy[i])
        atomT1.st()
        if slab == 1:
            atomsT_L.append(atomT1)
        if slab == 2:
            atomsT_R.append(atomT1)
        
    print(len(atomsT_L))
    t1.clear()
    t1.penup()
    t1.ht()
    t1.goto(-250,-260)
    t1.write("Left click to select an atom", align = "left", font = ("arial", 10))
    t1.goto(-250,-280)
    t1.write("Right click to unselect", align = "left", font = ("arial", 10))
    wn.update()
    wn.mainloop()
    
# create slab's bonds
def plot_bonds(slab):   
    wn.tracer(0)  
    if slab == 1:
        bondT1.clear()
        bondT1.color("blue")
        bond_status[0] = True
        
    if slab ==2:
        bondT2.clear()
        bondT2.color("red")
        bond_status[1] = True
    aR=1
    e=0.1
    neighbors = []

    try:
        if slab == 1:
            coordx = coordx_L
            coordy = coordy_L
            #coordz = coordz_L
        if slab == 2:
            coordx = coordx_R
            coordy = coordy_R
            #coordz = coordz_R

        for i in range(len(coordx)): #len(coordx)

            jlist = []
            for j in range(len(coordx)):
                d = math.sqrt((coordx[j]-coordx[i])**2+(coordy[j]-coordy[i])**2)
                if ((d>=(1-e)*aR) and (d<=(1+e)*aR)):
                    if i!=j:
                        jlist.append(j)

            neighbors.append(jlist)
    
        if slab == 1:
            coordpx = coordpx_L
            coordpy = coordpy_L
        if slab == 2:
            coordpx = coordpx_R
            coordpy = coordpy_R
        for i in range(len(neighbors)):
            for j in range(len(neighbors[i])):
                if slab == 1:
                    bondT1.penup()
                
                    bondT1.goto(coordpx[i],coordpy[i])
                    bondT1.pendown()
                    bondT1.goto(coordpx[neighbors[i][j]],coordpy[neighbors[i][j]])
                if slab == 2:
                    bondT2.penup()
                
                    bondT2.goto(coordpx[i],coordpy[i])
                    bondT2.pendown()
                    bondT2.goto(coordpx[neighbors[i][j]],coordpy[neighbors[i][j]])
       
        wn.update()                    

    except:
        print("Import slab first")
        
# Not  used  yet      
def drag_atom(t1):
    wn.listen()  
    t1.ondrag(t1.dragging,3)  #drag with the right key    
    wn.mainloop()
    
# displays  the rotation center  coordinates in the entry widget    
def display_rot_center():
        global x_click,y_click
        tk.Label(frame4, text = "Rotation Center").grid(row=3, column = 0,pady = 2.5)

        vx = tk.StringVar()
        vy = tk.StringVar()
        vx.set(x_click)
        vy.set(y_click)
        numx = tk.Entry(frame4,textvariable=vx,text = "%s" %(vx))
        numx.config(width=6)
        numx.grid(row=3, column=1)
        numy = tk.Entry(frame4,textvariable=vy,text = "%s" %(vy))
        numy.config(width=6)
        numy.grid(row=3, column=2)

# atom selection  by clicking mouse left button
def unselect(x,y):
    global slab_status, atomsT_L, atomsT_R,x_click,y_click
    #wn.setup(600,800,startx=680, starty=0)
    wn.tracer(0)
    x_click = 0.0
    y_click = 0.0
    display_rot_center()
    
    #print("x_c1=",x_click,"y_c1=",y_click)
    if slab_status[0]:
        for i in range(len(atomsT_L)):
            atomsT_L[i].fillcolor("blue")
            atomsT_L[i].st()
    if slab_status[1]:
        for i in range(len(atomsT_R)):
            atomsT_R[i].color("red")
            atomsT_R[i].fillcolor("")
            atomsT_R[i].st()
    wn.update()        
def select_atoms(x,y):

    global atomsT_L,atomsT_R
    global  slab_status,x_click,y_click
   
    #wn.setup(600,900,startx=680, starty=0)
    wn.tracer(0)
    x_click = x
    y_click = y
    #print("x_c2=",x_click,"y_c2=",y_click)
    display_rot_center()
    if slab_status[0]:
        size =atomsT_L[0].shapesize()
    if slab_status[1]:
        size =atomsT_R[0].shapesize()
    shift=size[0]* 10
    
    bool1 = False
    bool2 = False
    
    if slab_status[0]:
        for i in range(len(atomsT_L)):
        
            
            bool1 =(x >= coordpx_L[i] - shift and x <= coordpx_L[i] + shift) and (y >= coordpy_L[i] - shift and y <= coordpy_L[i] +shift)
            
            if bool1:
                
                atomsT_L[i].fillcolor("yellow")
                atomsT_L[i].st()
                
    if slab_status[1]:
        for i in range(len(atomsT_R)):
            bool2 = (x >= coordpx_R[i] - shift and x <= coordpx_R[i] + shift) and (y >= coordpy_R[i] - shift and y <= coordpy_R[i] + shift)
            
            if bool2:
                
                atomsT_R[i].fillcolor("red")
                atomsT_R[i].st()
    wn.update()        
    
            

# mouse binding for atom selection

def click_main():
    
    wn.onclick(select_atoms,1)
    wn.onclick(unselect,3,add = True)
    wn.mainloop()
  
# Calculate and display the coincidence lattice
def coincidence_lattice():

    global midpointsx,midpointsy
    aR=1
    e=15                           
    midpoints_px=[]
    midpoints_py=[]
    midpointsx=[]
    midpointsy=[]
    k = 0
    try:
        
        
        file1 = open("midpoints.txt", "w")
        for iL in range(len(coordx_L)):
            for iR in range(len(coordx_R)):
                
                d = math.sqrt((coordx_L[iL]-coordx_R[iR])**2+(coordy_L[iL]-coordy_R[iR])**2)
                
                if d < (aR*e/100):
##                    print
##                    print
##                    print("d_midpoint=",d)
                    midpointsx.append((coordx_L[iL]+coordx_R[iR])/2.0)
                    midpointsy.append((coordy_L[iL]+coordy_R[iR])/2.0)
                    midpoints_px.append(round((coordpx_L[iL]+coordpx_R[iR])/2.0,1))
                    midpoints_py.append(round((coordpy_L[iL]+coordpy_R[iR])/2.0,1))
                    k+=1
        #midpoints_px,midpoints_py = convert_angstrom_pixels(midpointsx,midpointsy,alpha,betax,betay)
        print("num coinc pt=",k)
        wn.tracer(0)
        if len(midpoints_px) > 0:
            
            for i in range(len(midpoints_px)):       
                atomC.goto(midpoints_px[i],midpoints_py[i])
                atomC.stamp()
                
        wn.update()
        for i in range(len(midpointsx)):
            
            rec = "%0.4f,%0.4f " %(midpointsx[i],midpointsy[i])
            file1.write(rec+"\n")
        file1.close()
    except:
        print(" upload slabs first")
# Clear coincidence lattice
def clear_coincidence_Lat():
    atomC.clearstamps()

#Help
def About():
    root = tk.Tk()
    root.title(" About Nanojunction Builder")
    print("""Nano Junction allows you to build the junction between
             two overlaping nanostructures.The current version built
             the junction between two sheets of graphene""")
    S = tk.Scrollbar(root)
    T = tk.Text(root, height=4, width=50)
    S.pack(side=tk.RIGHT, fill=tk.Y)
    T.pack(side=tk.LEFT, fill=tk.Y)
    S.config(command=T.yview)
    T.config(yscrollcommand=S.set)
    quote = """Nano Junction allows you to build the junction between two overlaping nanostructures. The current version built the junction between two sheets of graphene"""
    T.insert(tk.END, quote)

# Display translation menu    
def Translate_menu():
    
    symbols = ["\u21E7",'\u21E6',"\u21E8","\u21E9"]
    # initializing of the  choice1(slab1, slab2 or union
    v1 = tk.IntVar()
    v1.set(0)
    # initializing of  choice2, distance or screen
    v2= tk.IntVar()
    v2.set(0)  


    choices1 = [
        "Slab1 ",
        "Slab2",
        "Union"
    ]
    choices2 = [
        "Distance (Å) ",
        "Screen (%)"
    ]

    def ShowChoice_1():
        print(v1.get())

    i = 1

    for val,choice in enumerate(choices1):
        tk.Radiobutton(frame3, 
                      text=choice,
                      padx = 20, 
                      variable=v1, 
                      command=ShowChoice_1,
                      value =val).grid(row = i, column = 0,sticky = tk.W)  
        i += 1     
    
    def ShowChoice_2():
        print(v2.get())
    
    
    i = 1

    for val,choice in enumerate(choices2):
        
        tk.Radiobutton(frame3, 
                      text=choice,
                      padx = 20, 
                      variable=v2, 
                      command=ShowChoice_2,
                    value =val).grid(row = i, column = 1,sticky = tk.W)
        
        i += 1     

    

    item_1 = tk.IntVar()
    item_2 = tk.IntVar()
    item_1.set(0)
    item_2.set(0)
    def get_item_values(v):
        global dist,direction

        if v.get() == 0:
            dist = float(item_1.get())
            print(dist)
        if v.get() == 1:
            dist = float(item_2.get())
            print(dist)


    def left_button():
        global dist,direction
        direction = "left"
        print(direction)
        translate_atom(dist,direction,v1,v2)
    def right_button():
        global dist,direction
        direction = "right"
        print(direction)
        translate_atom(dist,direction,v1,v2)
    def up_button():
        global dist,direction
        direction = "up"
        print(direction)
        translate_atom(dist,direction,v1,v2)
    def down_button():
        global dist,direction
        direction = "down"
        print(direction)
        translate_atom(dist,direction,v1,v2)
    


    btn1 = tk.Button(frame3, text=symbols[0],fg="black",background = "light grey",command =up_button )
    btn1.grid(row = 1, column = 4)    

    btn2 = tk.Button(frame3, text=symbols[1],fg="black",background = "light grey",command =left_button )
    btn2.grid(row = 2,column =3,padx = 1)

    btn3 = tk.Button(frame3, text=symbols[2],fg="black",background = "light grey", command =right_button)
    btn3.grid(row = 2,column = 5,padx =1)
    btn4 = tk.Button(frame3, text=symbols[3],fg="black",background = "light grey", command =down_button)
    btn4.grid(row = 3, column = 4)
    
    item_1 = tk.Spinbox(frame3, from_= 0, to = 250, increment = 0.1, width = 5,command = lambda:get_item_values(v2))
    item_1.grid(row = 1, column = 2,padx =5)
    item_2 = tk.Spinbox(frame3, from_= 0, to = 100, increment = 0.5, width = 5,command = lambda:get_item_values(v2))
    item_2.grid(row = 2, column = 2,padx =5)
     
# Translation of the sheet(s)
def translate_atom(dist,direction,v1,v2):

    global coordx_L,coordy_L,coordx_R,coordy_R,coordpx_L,coordpy_L,coordpx_R,coordpy_R
    
    global beta_L,beta_R,alpha
    global x1,x2,xwin1,xwin2
    global y1,y2,ywin1,ywin2
    global xPix1,xPix2,yPix1, yPix2
    global atomsT_L,atomsT_R,slab_status,bond_status
    
    wn.tracer(0)

    coordpx_trans = []
    coordpy_trans = []
    
    
    
    coordx_trans =[]
    coordy_trans  = []

    
    if direction == "left":
        dist = -dist
    if direction == "down":
        dist = -dist
    
    print("dist =", dist)

    
    # Translate slab1  if choice = 0 (slab1) or 2 (union)         
    
    if v1.get() == 0 or v1.get() == 2:
        if slab_status[0]:
            coordpx_trans = []
            coordpy_trans = []
            
            coordx_trans =[]
            coordy_trans  = []
            xPix1 = xwin1 + beta_L[0] *(xwin2-xwin1)  
            yPix1 = ywin1 + beta_L[1] *(ywin2-ywin1)


            xPix2 = xPix1 + alpha *(x2-x1)
            yPix2 = yPix1 + alpha *(y2-y1)
            # translate by angstrom
            if v2.get() == 0:
                dist_pix = (xPix2-xPix1)*dist/(x2-x1)
                dist_piy = (yPix2-yPix1)*dist/(y2-y1)
                # Translate Left slab if choice = 0 or 2
                for i in range(len(coordx_L)):
                    if direction == "left" or direction == "right":
                        coordx_trans.append(coordx_L[i] + dist)
                        coordpx_trans.append(coordpx_L[i] + dist_pix)
                        coordy_trans.append(coordy_L[i])
                        coordpy_trans.append(coordpy_L[i])
                    elif direction == "up" or direction == "down":
                        coordy_trans.append(coordy_L[i] + dist)
                        coordpy_trans.append(coordpy_L[i] + dist_piy)
                        coordx_trans.append(coordx_L[i])
                        coordpx_trans.append(coordpx_L[i])
                
            #translate by a percentage
            elif v2.get() == 1:
                dist_angstrom_x = (dist/100)*(x2-x1)
                dist_angstrom_y = (dist/100)*(y2-y1)
                dist_pix = (xPix2-xPix1)*dist_angstrom_x/(x2-x1)
                dist_piy = (yPix2-yPix1)*dist_angstrom_y/(y2-y1)
                for i in range(len(coordx_L)):
                    if direction == "left" or direction == "right":
                        coordx_trans.append(coordx_L[i] + dist_angstrom_x)
                        coordpx_trans.append(coordpx_L[i] + dist_pix)
                        coordy_trans.append(coordy_L[i])
                        coordpy_trans.append(coordpy_L[i])

                    elif direction == "up" or direction == "down":
                        coordy_trans.append(coordy_L[i] + dist_angstrom_y)
                        coordpy_trans.append(coordpy_L[i] + dist_piy)
                        coordx_trans.append(coordx_L[i])
                        coordpx_trans.append(coordpx_L[i])
            if direction == "left" or direction == "right":            
                beta_L=[beta_L[0]+ dist_pix/(xwin2-xwin1),beta_L[1]]
            if direction == "up" or direction == "down":
                beta_L=[beta_L[0],beta_L[1]+ dist_piy/(ywin2-ywin1)]
            print( "beta Left =", beta_L)
            coordx_L = coordx_trans
            coordpx_L = coordpx_trans
            coordy_L = coordy_trans
            coordpy_L = coordpy_trans
           

            #print("1st slab",len(atomsT_L))
            for i in range(len(atomsT_L)):
                atomsT_L[i].clear()
                atomsT_L[i].goto(coordpx_trans[i],coordpy_trans[i])
                atomsT_L[i].st()
             
            
            if bond_status[0]:
              
                bondT1.clear()
                plot_bonds(1)
                
            wn.update()
            
    # Translate slab2  if choice = 1 (slab2) or 2 (union)        
    if v1.get() == 1 or v1.get() == 2:
        if slab_status[1]:
            coordpx_trans = []
            coordpy_trans = []
            
            coordx_trans =[]
            coordy_trans  = []
            xPix1 = xwin1 + beta_R[0] *(xwin2-xwin1)  
            yPix1 = ywin1 + beta_R[1] *(ywin2-ywin1)


            xPix2 = xPix1 + alpha *(x2-x1)
            yPix2 = yPix1 + alpha *(y2-y1)
            
            # translate by angstrom
            if v2.get() == 0:
                
                dist_pix = (xPix2-xPix1)*dist/(x2-x1)
                dist_piy = (yPix2-yPix1)*dist/(y2-y1)
                # Translate Left slab if choice = 0 or 2
                for i in range(len(coordx_R)):
                    if direction == "left" or direction == "right":
                        coordx_trans.append(coordx_R[i] + dist)
                        coordpx_trans.append(coordpx_R[i] + dist_pix)
                        coordy_trans.append(coordy_R[i])
                        coordpy_trans.append(coordpy_R[i])
                    elif direction == "up" or direction == "down":
                        coordy_trans.append(coordy_R[i] + dist)
                        coordpy_trans.append(coordpy_R[i] + dist_piy)
                        coordx_trans.append(coordx_R[i])
                        coordpx_trans.append(coordpx_R[i])
                
               
            #translate by a percentage
            elif v2.get() == 1:
                dist_angstrom_x = (dist/100)*(x2-x1)
                dist_angstrom_y = (dist/100)*(y2-y1)

                dist_pix = (xPix2-xPix1)*dist_angstrom_x/(x2-x1)
                dist_piy = (yPix2-yPix1)*dist_angstrom_y/(y2-y1)
                
                
                #print("d1=",dist_pix,"d2=",dist_piy)

                for i in range(len(coordx_R)):
                    if direction == "left" or direction == "right":
                        coordx_trans.append(coordx_R[i] + dist_angstrom_x)
                        coordpx_trans.append(coordpx_R[i] + dist_pix)
                        coordy_trans.append(coordy_R[i])
                        coordpy_trans.append(coordpy_R[i])

                    elif direction == "up" or direction == "down":
                        coordy_trans.append(coordy_R[i] + dist_angstrom_y)
                        coordpy_trans.append(coordpy_R[i] + dist_piy)
                        coordx_trans.append(coordx_R[i])
                        coordpx_trans.append(coordpx_R[i])
            if direction == "left" or direction == "right":            
                beta_R=[beta_R[0]+ dist_pix/(xPix2-xPix1),beta_R[1]]
            if direction == "up" or direction == "down":
                beta_R=[beta_R[0],beta_R[1]+ dist_piy/(xPix2-xPix1)]
            print( "beta Right =", beta_R)            
            coordx_R = coordx_trans
            coordpx_R = coordpx_trans
            coordy_R = coordy_trans
            coordpy_R = coordpy_trans
       

            for i in range(len(atomsT_R)):
                atomsT_R[i].clear()
                atomsT_R[i].goto(coordpx_trans[i],coordpy_trans[i])
                atomsT_R[i].st()

            
            if bond_status[1]:
                bondT2.clear()                
                plot_bonds(2) 
            wn.update()

#displays rotation menu            
def Rotation_menu():
    symbols = ['\u21B7','\u21B6']

    slab= tk.IntVar()
    slab.set(0)  # initializing the choice1 slab


    choices1 = [
        "Slab1 ",
        "Slab2"
    ]
    
    def ShowChoice_1():
        print(slab.get())

    def get_angle_value():
        global angle
        angle = float(item_1.get())
        print(angle)
        
    
    i = 1

    for val,choice in enumerate(choices1):
        
        tk.Radiobutton(frame4, 
                      text=choice,
                      padx = 20, 
                      variable=slab, 
                      command=ShowChoice_1,
                    value =val).grid(row = i, column = 0,sticky = tk.W)
        
        i += 1
    
    rot_dir = " "    
    def clock_wise_button():
        global angle,rot_dir
        rot_dir = "clockw"
        angle = float(item_1.get())
        print(rot_dir,angle)
        rotate_slab(angle,slab)
        
    def couterclock():
        global angle, rot_dir
        rot_dir = "aclockw"
        angle=-float(item_1.get())
        print(rot_dir,angle)
        rotate_slab(angle,slab)
        
    btn1 = tk.Button(frame4, text=symbols[0],fg="black",background = "light grey",command =clock_wise_button)
    btn1.grid(row = 1, column = 2)    

    btn2 = tk.Button(frame4, text=symbols[1],fg="black",background = "light grey",command =couterclock)
    btn2.grid(row = 2,column =2,padx = 1)
    item_1 = tk.Spinbox(frame4, from_= 0, to = 180, increment = 0.1, width = 5,command = get_angle_value)
    item_1.grid(row = 1, column =1,padx =5)
    #cx_px,cy_px =show_selected_rot_center()
    display_rot_center()
    
#rotate a slab about an atom    


def rotate_slab(angle,slab):
    global coordx_L,coordy_L,coordx_R,coordy_R,coordpx_L,coordpy_L,coordpx_R,coordpy_R
    global rot_dir
    global beta_L,beta_R,alpha
    global x1,x2,xwin1,xwin2
    global y1,y2,ywin1,ywin2
    global xPix1,xPix2,yPix1, yPix2
    global atomsT_L,atomsT_R,slab_status,x_click,y_click
   
    #get center of rotation coordinates
    
    cx_px,cy_px = x_click,y_click
    
    #print("cx_px=", cx_px, "cy_px=", cy_px)
    #print("beta_L=", beta_L, "beta_R=", beta_R)
    
    cx_px = float(cx_px)
    cy_px = float(cy_px)

    # Rotate Slab1
    if slab.get() == 0 and slab_status[0]:         
         coordpx_rot = []
         coordpy_rot = []
         coordx_rot =[]
         coordy_rot= []
         print(angle, "angle")
         angle_rad = angle*math.pi/180

         cx = (cx_px - xPix1)*(x2-x1)/(xPix2-xPix1)+x1
         cy = (cy_px - yPix1)*(y2-y1)/(yPix2-yPix1)+y1
         print("cc", cx,cy)
         for i in range(len(coordx_L)):
            
            xPix = cx_px+(coordpx_L[i]-cx_px)*math.cos(angle_rad) +(coordpy_L[i]-cy_px)*math.sin(angle_rad)
            x=cx+(coordx_L[i]-cx)*math.cos(angle_rad) +(coordy_L[i]-cy)*math.sin(angle_rad)
            coordx_rot.append(x)
            coordpx_rot.append(xPix)
            y = cy-(coordx_L[i]-cx)*math.sin(angle_rad) +(coordy_L[i]-cy)*math.cos(angle_rad)
            yPix = cy_px-(coordpx_L[i]-cx_px)*math.sin(angle_rad) +(coordpy_L[i]-cy_px)*math.cos(angle_rad)
            coordy_rot.append(y)
            coordpy_rot.append(yPix)
            
         #coordpx_rot,coordpy_rot = convert_angstrom_pixels(coordx_rot,coordy_rot,alpha,beta_L)
         coordx_L = coordx_rot
         coordy_L = coordy_rot
         coordpx_L = coordpx_rot
         coordpy_L = coordpy_rot

         
         wn.tracer(0)

         #print("rotation",len(atomsT_L))
         for i in range(len(atomsT_L)):
            atomsT_L[i].clear()
            atomsT_L[i].goto(coordpx_rot[i],coordpy_rot[i])
            atomsT_L[i].st()
         if bond_status[0]:    
            bondT1.clear()
            plot_bonds(1)             
         wn.update()

    # Rotate Slab2
    if slab.get() == 1 and slab_status[1]:         
         coordpx_rot = []
         coordpy_rot = []
         coordx_rot =[]
         coordy_rot= []
         #print(angle, "angle")
         
         # Convert angle into Radiant
         angle_rad = angle*math.pi/180
         xPix1 = xwin1 + beta_R[0] *(xwin2-xwin1)  
         yPix1 = ywin1 + beta_R[1] *(ywin2-ywin1)
         xPix2 = xPix1 + alpha *(x2-x1)
         yPix2 = yPix1 + alpha *(y2-y1)
         cx = (cx_px - xPix1)*(x2-x1)/(xPix2-xPix1)+x1
         cy = (cy_px - yPix1)*(y2-y1)/(yPix2-yPix1)+y1
         #print("cc", cx,cy)
         for i in range(len(coordx_R)):
            
            xPix = cx_px+(coordpx_R[i]-cx_px)*math.cos(angle_rad) +(coordpy_R[i]-cy_px)*math.sin(angle_rad)
            x=cx+(coordx_R[i]-cx)*math.cos(angle_rad) +(coordy_R[i]-cy)*math.sin(angle_rad)
            coordx_rot.append(x)
            coordpx_rot.append(xPix)
            y = cy-(coordx_R[i]-cx)*math.sin(angle_rad) +(coordy_R[i]-cy)*math.cos(angle_rad)
            yPix = cy_px-(coordpx_R[i]-cx_px)*math.sin(angle_rad) +(coordpy_R[i]-cy_px)*math.cos(angle_rad)
            coordy_rot.append(y)
            coordpy_rot.append(yPix)
            

         coordx_R = coordx_rot
         coordy_R = coordy_rot
         coordpx_R = coordpx_rot
         coordpy_R = coordpy_rot

         
         wn.tracer(0)

         #print("rotation",len(atomsT_R))
         for i in range(len(atomsT_R)):
            atomsT_R[i].clear()
            atomsT_R[i].goto(coordpx_rot[i],coordpy_rot[i])
            atomsT_R[i].st()
         if bond_status[0]:
              
             bondT2.clear()
             plot_bonds(2)
                

         wn.update()

 

# Quit ; Clear screen
def _quit(): 
    window.quit()     # stops mainloop
    window.destroy()# this is necessary on Windows to prevent
    wn.bye()                # Fatal Python Error: PyEval_RestoreThread: NULL tstate

def clearscreen():
    global num_slab,slab_status,coordx_L,coordy_L,coordx_R,coordy_R,coordz_R,coordpx_L
    global coordpy_L,coordpz_L, coordpx_R,coordpy_R,coordpz_R
    global atoms_sel,atomsT_L,atomsT_R,atoms_L,atoms_R,atomsT, new_Lattice
    global midpoints_px,midpoints_py,atomT,reset,x1_y1_flag,beta_L,beta_R,bond_status
   

    slab_status = [False]*num_slab
    bond_status = [False]*num_slab
    midpoints_px =[]
    midpoints_py=[]
    coordx_L = []
    coordy_L = []
    coordz_L = []
    coordx_R = []
    coordy_R = []
    coordz_R = []
    coordpx_L = []
    coordpy_L = []
    coordpz_L = []

    coordpx_R = []
    coordpy_R = []
    coordpz_R = []
    coordx_L_trans = []
    coordy_L_trans = []
    coordz_L_trans = []
    coordx_R_trans = []
    coordy_R_trans = []
    coordz_R_trans = []
    atoms_sel = []
    atomsT_L = []
    atomsT_R=[]
    atoms_L = []
    atoms_R = []
    atomsT =[]
    new_Lattice = []
    direction = " "
    x1_y1_flag = False
    beta_L=[0,0.2]
    beta_R=[0,0.2]
    wn.clearscreen()
    Translate_menu()
    Rotation_menu()
    bondT1.clear()
    bondT2.clear()
    click_main()
   
   
# Defining Frames, Canvas, Buttons and binding Keys and events    

frame1 = tk.Frame(master=window, width=600, height=350, bg="light grey",relief = tk.SUNKEN,borderwidth=2)



frame1.place(x=10, y=5)
canvas1 = tk.Canvas(master = frame1, width = 400, height = 350)

canvas1.pack()


frames = []
buttons = []
for j in range(3):
   frame =  tk.Frame(master=window, width=248,height=157, bg="light grey",relief = tk.SUNKEN,borderwidth=2)
 
   
   frame.place(x=420, y=5+j*160)
   frames.append(frame)
frame2 = tk.Frame(master=window, width=400,height=200, bg="light grey",relief = tk.SUNKEN,borderwidth=1)
frame2.pack(fill=tk.BOTH, side=tk.LEFT)
frame2.place(x=10, y=360)
step = 4
for i in range(4):
    for j in range(4):
        if i == 0 and j == 0:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = lambda: plot_sheets(alpha,beta_L,1))
        elif i == 0 and j == 1:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = lambda: plot_sheets(alpha,beta_R,2))
        elif  i == 1 and j == 0:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = lambda: plot_bonds(1)) 
        elif  i == 1 and j == 1:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = lambda: plot_bonds(2))
        elif  i == 1 and j == 2:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = coincidence_lattice)
        elif  i == 1 and j == 3:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = clear_coincidence_Lat)    
            
        elif  i == 0 and j == 3:
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey", command = export_file)
        else:    
            btn = tk.Button(frame2, text=texts[i*step+j],width=12,fg="black",activebackground = "dark grey")  
        btn.grid(row=i, column=j, padx=3, pady=2)
        buttons.append(btn)
    
    

frame3 = tk.Frame(master=window, width=550,height=200,relief = tk.SUNKEN,borderwidth=2)

frame3.place(x=10, y=500)
frame4 = tk.Frame(master=window, width=270,height=200,relief = tk.SUNKEN,borderwidth=2)

frame4.place(x=380, y=500)


label1 = tk.Label(master=frame3, anchor = tk.CENTER,text="Translation", font =("arial",11,"bold"))


label1.grid(sticky =tk.NW, row = 0,column = 0)
label2 = tk.Label(master=frame4, anchor = tk.CENTER,text="Rotation", font =("arial",11,"bold"))


label2.grid(sticky =tk.NW, row = 0,column = 0)

fig = Figure(figsize=(2.14, 1.34), dpi=115)
t = np.arange(0, 2, .01)
fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
canvas = FigureCanvasTkAgg(fig, master=frames[0])  
canvas.draw()
canvas.get_tk_widget().pack() #fill=tk.BOTH, expand=1

#toolbar = NavigationToolbar2Tk(canvas, frames[0])

##toolbar.SetToolBitmapSize(12,12)

#toolbar.update()
fig = Figure(figsize=(4, 3), dpi=100)
t = np.arange(0, 3, .01)
fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
canvas2 = FigureCanvasTkAgg(fig, master=canvas1)  # A tk.DrawingArea.
canvas2.draw()
canvas2.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

toolbar = NavigationToolbar2Tk(canvas2, canvas1)
toolbar.update()
#print(toolbar.toolitems) 

    
menu = tk.Menu(window)
window.config(menu=menu)

helpmenu = tk.Menu(menu)
menu.add_cascade(label="Menu", menu=helpmenu)
helpmenu.add_command(label="About...", command=About)
helpmenu.add_separator()
helpmenu.add_command(label="Topics", command=helpmenu.quit)

helpmenu.add_separator()
helpmenu.add_command(label="Clear screen", command=clearscreen)
helpmenu.add_separator()
helpmenu.add_command(label="Quit", command=_quit)


if __name__ == "__main__":
    Translate_menu()
    Rotation_menu()
    click_main()
    window.mainloop()


































