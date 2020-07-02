# %matplotlib notebook
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits
from mpl_toolkits import mplot3d
import re
import matplotlib.animation

#time = []
time = np.linspace(0,20,1000)
pos=[]
pos_temp =[]
f = open("pos.txt","r")
i=0
for lines in f:
    if lines.startswith("t"):
        tline = lines.strip()
        t=float(''.join(re.findall('\d*\.?\d+',tline)))
        #time.append(t)
    elif lines.startswith("\n"):
       #print(pos_temp)
        pos.append(pos_temp)
        pos_temp =[]

    else:
        lines=lines.strip()
        lines= lines.split()
        lines = list(map(float, lines))
        pos_temp.append(lines)

def animation(p):
    current_X = []
    current_Y = []
    current_Z = []
    current_time = time[p]
    current_positions = pos[p]
    for q in range(0,len(current_positions)):
        current_X.append(pos[p][q][0])
        current_Y.append(pos[p][q][1])
        current_Z.append(pos[p][q][2])

    if (p!=0):
        graph._offsets3d = (current_X,current_Y,current_Z)
    else:
        pass

    return current_X,current_Y,current_Z
    print(current_time,current_positions,"\n",current_X,current_Y,current_Z,"\n")


fig =plt.figure()
ax =fig.add_subplot(111,projection='3d')
p = np.arange(1,len(time))
init_graph = animation(0)
graph = ax.scatter(init_graph[0],init_graph[1],init_graph[2],c='r',marker='o')


ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
ax.set_zlabel('z axis')
ax.set_xlim(-100, 100)
ax.set_ylim(-100, 100)
ax.set_zlim(-100, 100)
ani = matplotlib.animation.FuncAnimation(fig,animation,p,interval=10, blit=True)
plt.show()
