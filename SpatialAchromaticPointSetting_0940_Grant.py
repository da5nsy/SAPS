from psychopy import visual, event, core, data 
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

#set up definitions
m,n = [1366,768] #define screen size (tablet =[1366,768] laptop =[1920,1080])
MdpX=np.int16(m/6) #max displacment X
MdpY=np.int16(n/6) #max displacment Y
k=np.sqrt(((m/2)**2)+((m/2)**2)) 
stimSize=np.int16(m+(2*MdpY)+(2*(k-(m/2))))+1
timeAllowed=100 #sets time allowed for input
sampleInterval=0.1 #sets accuracy of timing 

basicStim = np.array([[ 1, -0.5],[ -0.5, 1],])

theWin   = visual.Window(screen=0, 
    size=[m,n], 
    color=(-1,-1,-1),
    monitor='testMonitor',
    fullscr=
    True,
    #False,
    units='pix')

theStim = visual.ImageStim(theWin,
    image=
    #b,
    "zazzle_60_50_8bit.tif",
    #basicStim, 
    #theGrad, 
    pos=(0,0), 
    size=
    #(400,400),
    (stimSize,stimSize),
    units='pix', 
    interpolate=False,
    autoLog=False,)

#boring gubbins
theMouse=event.Mouse(visible=False)#defines mouse, and makes invisible
originalPos=theMouse.getPos() #needs to be called so that mouseMoved has something to compare against
timeElapsed=0 #setup value for while loop

variations=200
checker_switch=variations-11
repeats=[2,7]
dpX=np.random.random_integers(-MdpX, MdpX, variations)
#dpX=[99]*variations
dpX[repeats]=dpX[repeats[1]]

dpY=np.random.random_integers(-MdpY, MdpY, variations)
#dpY=[99]*variations
dpY[repeats]=dpY[repeats[1]]

ori=np.random.random_integers(0, 360, variations)
#ori=[0]*variations
ori[repeats]=ori[repeats[1]]

stimList=[]
for i in range(variations):
    stimList.append( {'dpX':dpX[i], 'dpY':dpY[i] ,'ori':ori[i]} )

#Data storage Pre
trials = data.TrialHandler(stimList,1, method='sequential'
,extraInfo= {'participant':"DG Basement 3",'hand':"R",'repeats':str(repeats),'checker_switch':checker_switch}
)
trials.data.addDataType('x-co')
trials.data.addDataType('y-co')
#...do I need to add other dataTypes? Seems to work fine without currently...
magnitude=np.zeros(variations)
direction=np.zeros(variations)
cnt=-1

#stumuli are go
#event.waitKeys() #waits for keypress
for thisTrial in trials:
    theStim.pos=(thisTrial['dpX'],thisTrial['dpY'])
    theStim.ori=(thisTrial['ori'])
    theStim.draw()
    theWin.flip()
    timeElapsed=0 #setup value for while loop
    cnt=cnt+1
    if cnt==checker_switch:
        theStim.image=basicStim 
        theStim.size=(300,300)
    if event.getKeys(keyList=['escape','q']):
        theWin.close()
        core.quit()
    while timeElapsed<timeAllowed:
        core.wait(sampleInterval) 
        timeElapsed=timeElapsed+sampleInterval
        if theMouse.mouseMoved(): 
            newPos=theMouse.getPos() #record position
            trials.data.add('x-co', newPos[0]) #add the data to our set
            trials.data.add('y-co', newPos[1]) 
            xCoCor=newPos[0]-thisTrial['dpX']
            trials.data.add('x-co-cor',xCoCor) 
            yCoCor=newPos[1]-thisTrial['dpY']
            trials.data.add('y-co-cor',yCoCor)
            magnitude[cnt]=np.sqrt((xCoCor**2)+(yCoCor**2))
            trials.data.add('magnitude',magnitude[cnt])#magnitude from (non-displaced) centre
            direction[cnt]=(np.arctan2(yCoCor,xCoCor)* 180 / np.pi)+thisTrial['ori']
            if direction[cnt]<0:
                direction[cnt]=360+direction[cnt]
            trials.data.add('direction',direction[cnt]# if direction[cnt]>0 else 360+direction[cnt])
            )
            trials.data.add('toc',timeElapsed)
            timeElapsed=timeAllowed #ends while loop, comment out if multiple selections to be allowed
theWin.close()

#Results

#resultsIm = np.zeros(variations)
#for i in range(0,variations):
#    resultsIm[i]=theGrad[stimSize/2,stimSize/2]
#
#print resultsIm 

#trials.printAsText(stimOut=['dpX','dpY','ori'], #write summary data to screen 
#                  dataOut=['x-co_raw','y-co_raw','magnitude_raw','direction_raw','toc_raw'])
#
trials.saveAsExcel(fileName='200s', 
                  sheetName = 'Data_'+str(core.getAbsTime()),
                  stimOut=['dpX','dpY','ori'], 
                  dataOut=['x-co_raw','y-co_raw','magnitude_raw','direction_raw','toc_raw']
                  )
                  

#Polar plot
r = magnitude
theta = direction / (180 /np.pi)
ax = plt.subplot(111, projection='polar')
ax.plot(theta, r, color='r', linewidth=1)
ax.set_rmax(max(r+50))
ax.grid(True)
#ax.set_title("A line plot on a polar axis", va='bottom')
plt.show()
