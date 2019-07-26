from psychopy import visual, event
import numpy as np
from linz import linz
from XYZtoRGB import XYZtoRGB

basicStim = np.array([
        [ 1, 1, 1, 1],
        [ 1, 1, 1, 1],
        [ 1, 1, 1, 1],
        [ 1, 1, 1, 1]
        ])

theWin   = visual.Window(screen=0, 
    color=(-1,-1,1),
    size=[1366/2,768],
    monitor='testMonitor',
    fullscr=False,
    units='pix',
    pos=[0,0]
    )

theStim = visual.ImageStim(theWin,
    image=basicStim,
    color=[-1,-1,-1], 
    size=(theWin.size))
    
theText=visual.TextStim(theWin,
    text='start',
    pos=(290, 50),
    ori=270)
    
#theStim.draw()
theText.draw()
theWin.flip()
event.waitKeys()

#working for extremes, but not for intermediaries, WTF?!
#Oh, because the XYZ values cannot be independently varied? How does that work...
#It's happy so long as the differences are not too great between the X,Y, and Z
#that's some cray s***
#theStim.color= (((np.array(linz(XYZtoRGB(270,270,270))))*2)-1)
#theStim.color=[-1,-1,1]
#theStim.draw()
#theWin.flip()
#event.waitKeys()

#print theStim.color
#print (np.array(linz(XYZtoRGB(350,380,430)))*2)-1
#print (np.array(linz(XYZtoRGB(0,0,0)))*2)-1


for i in range(0,256,15):
    norm=(((float(i))/255)*2)-1
    theStim.color=[norm,-1,-1]
    theText.text=str(i)
    theStim.draw()
    theText.draw()
    theWin.flip()
    event.waitKeys()
#
#linearized
#for i in range(0,256,15):
#    norm=(float(i))/255
#    norm = ((linz1(norm))*2)-1
#    theStim.color=[-1,-1,norm]
#    theText.text=str(i)
#    theStim.draw()
#    theText.draw()
#    theWin.flip()
#    event.waitKeys()

#XYZ + linearized
#for i in range(0,256,15):
#    norm=(float(i))/255
#    norm = ((linz1(norm))*2)-1
#    theStim.color=(((np.array(linz(XYZtoRGB(i,i,i))))*2)-1)
#    theText.text=str(i)
#    theStim.draw()
#    theText.draw()
#    theWin.flip()
#    event.waitKeys()
#
