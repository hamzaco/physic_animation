import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.transforms as tr
import random
class Particle:
    def __init__(self,center=(0,0),radius=1,velocity=np.array([0,0])):
        self.object=patches.Circle(center,radius,color='k')
        self.position=center
        self.velocity=velocity
        self.mass=radius**2
    def get_radius(self):
        return self.object.get_radius()

    def get_center(self):
        return self.object.get_center()
    def set_center(self,center):
        self.object.set_center(center)

    def get_boundaries(self):
        angles=np.arange(0,2*np.pi,0.001)
        boundaries=[]
        r=self.get_radius()
        c=self.get_center()
        for angle in angles:
            boundaries.append((c[0]+r*np.sin(angle),c[1]+np.cos(angle)))
        return boundaries
    def is_inside(self,box):
        xlim=box.intervalx
        ylim=box.intervaly
        tempx= self.get_center()[0]+self.get_radius()
        tempy= self.get_center()[1]+self.get_radius()
        if xlim[0] > self.get_center()[0]-self.get_radius():
            return False,np.array([1,0])
        if xlim[1] < tempx:
            return False,np.array([-1,0])
        if ylim[0] > self.get_center()[1]-self.get_radius():
            return False,np.array([0,1])
        if ylim[1] < tempy:
            return False,np.array([0,-1])
        return True, None

    def bounce(self,normal):
        if normal[0]==0:
            self.velocity[1]=-self.velocity[1]
        else:
            self.velocity[0]=-self.velocity[0]

    def add(self,ax):
        ax.add_patch(self.object)

    def move(self,dt,Box,force='magnetic'):
        c=self.get_center()
        temp=self.is_inside(Box)
        if not temp[0]:
            self.bounce(temp[1])

        self.object.set_center(c+self.velocity*dt)
        if force=='magnetic':
            force=np.array([self.velocity[1],-self.velocity[0]])/np.linalg.norm(self.velocity)*0.01
            self.velocity = self.velocity - force / self.mass * dt

    def get_direction(self, particles, r0):
        neighbours = []
        for op in particles:
            if r0 > tuple_distance(op.position, self.position):
                neighbours.append(op.velocity)


def length(a):
    return np.sqrt(np.dot(a,a))

def angle_between(a,b):
    return np.arccos(np.dot(a,b)/(length(a)*length(b)))

def tuple_distance(t1,t2):
    v=np.array(t1)-np.array(t2)
    return np.sqrt(v[0]**2+v[1]**2)







def random_particles(n,ax):
    particles=[]
    for i in range(n):
        p=Particle(center=(random.randrange(10,40,1)/10,random.randrange(10,40,1)/10),
                   velocity=np.array([random.random()*10,random.random()*10]),
                   radius=random.randrange(20,40,1)/10000)
        particles.append(p)
        ax.add_patch(p.object)
    return particles
def move_particles(particles,dt,box,force):
    for p in particles:
        p.move(dt,box,force)

fig, ax = plt.subplots()

t = np.arange(0.0, 5, 0.001)
plt.axis("off")
ax.set_xlim([0,5])
ax.set_ylim([0,5])
xlim = ax.get_xlim()
ylim = ax.get_ylim()
ax.set_aspect('equal', adjustable='datalim')
Box=plt.Polygon([(0,0),(0,5),(5,5),(5,0)],fill=False)
ax.add_patch(Box)
Box=tr.Bbox([[0,0],[5,5]])
particles=random_particles(100,ax)

def animate(i):
    dt=0.001
    move_particles(particles,dt,Box,'magnetic')
    return list(map(lambda x: x.object,particles))

# create animation using the animate() function
myAnimation = animation.FuncAnimation(fig, animate, \
                                      frames=1000, blit=True, repeat=True)

myAnimation.save('sine.html', fps=30, extra_args=['-vcodec', 'libx264'])