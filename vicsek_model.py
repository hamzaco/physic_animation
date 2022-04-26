import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.transforms as tr
import random


class Particle:
    def __init__(self, center=np.array([0, 0]), radius=1, orientation=np.array([1,0])):
        self.object = patches.Circle(center, radius, color='k')
        self.position = center
        self.orientation = orientation
        self.mass = radius ** 2

    def get_radius(self):
        return self.object.get_radius()

    def get_center(self):
        return self.object.get_center()

    def set_center(self, center):
        self.object.set_center(center)

    def is_inside(self, box):
        xlim = box.intervalx
        ylim = box.intervaly
        tempx = self.get_center()[0] + self.get_radius()
        tempy = self.get_center()[1] + self.get_radius()
        if xlim[0] > self.get_center()[0] - self.get_radius():
            return False, np.array([1, 0])
        if xlim[1] < tempx:
            return False, np.array([-1, 0])
        if ylim[0] > self.get_center()[1] - self.get_radius():
            return False, np.array([0, 1])
        if ylim[1] < tempy:
            return False, np.array([0, -1])
        return True, None

    def bounce(self, normal):
        c=self.object.get_center()
        o=self.orientation
        if normal[0] == 0 and normal[1]==1:
            nc=np.array([c[0],ylim[1]])
            self.orientation=np.array([o[0],abs(o[1])])
            #self.object.set_center(nc)
        elif normal[0] == 0 and normal[1]==-1:
            nc=np.array([c[0],ylim[0]])
            self.orientation=np.array([o[0],-abs(o[1])])
            #self.object.set_center(nc)
        elif normal[0] == 1 and normal[1]==0:
            nc=np.array([xlim[1],c[1]])
            self.orientation=np.array([abs(o[0]),o[1]])
            #self.object.set_center(nc)
        else:
            nc=np.array([xlim[0],c[1]])
            self.orientation=np.array([-abs(o[0]),o[1]])
            #self.object.set_center(nc)
    def add(self, ax):
        ax.add_patch(self.object)

    def get_direction(self, particles, r0):
        orientation_sum = np.array([0,0])
        for op in particles:
            if r0 > tuple_distance(op.position, self.position) > 0:
                orientation_sum = orientation_sum+op.orientation
        if sum(orientation_sum == 0)!=2:
            orientation_sum = orientation_sum/np.linalg.norm(orientation_sum)
            orientation_sum += np.array([random.random(), random.random()])/3
            orientation_sum = orientation_sum/np.linalg.norm(orientation_sum)
            return orientation_sum/np.linalg.norm(orientation_sum)
        return self.orientation
    def move(self, dt, v0, r0, Box,particles):
        temp = self.is_inside(Box)
        if not temp[0]:
            self.bounce(temp[1])
        c = self.get_center()
        self.object.set_center(c +dt*v0*self.orientation)

def length(a):
    return np.sqrt(np.dot(a, a))


def angle_between(a, b):
    return np.arccos(np.dot(a, b) / (length(a) * length(b)))


def tuple_distance(t1, t2):
    v = np.array(t1) - np.array(t2)
    return np.sqrt(v[0] ** 2 + v[1] ** 2)


def random_particles(n, ax):
    particles = []
    for i in range(n):
        o=np.array([random.random()-0.5, random.random()-0.5])
        p = Particle(center=np.array((random.randrange(10, 40, 1) / 10, random.randrange(10, 40, 1) / 10)),
                     orientation=o/np.linalg.norm(o),
                     radius=random.randrange(20, 40, 1) / 10000)
        particles.append(p)
        ax.add_patch(p.object)
    return particles


def move_particles(particles, dt, v0, r0, box,):
    directions=[]
    for p in particles:
        directions.append(p.get_direction(particles,r0))
    for i in range(0,len(particles)):
        particles[i].orientation= directions[i]
    for p in particles:
        p.move(dt, v0, r0, box, particles)


fig, ax = plt.subplots()

t = np.arange(0.0, 5, 0.001)
plt.axis("off")
ax.set_xlim([0, 5])
ax.set_ylim([0, 5])
xlim = ax.get_xlim()
ylim = ax.get_ylim()
ax.set_aspect('equal', adjustable='datalim')
Box = plt.Polygon([(0, 0), (0, 5), (5, 5), (5, 0)], fill=False)
ax.add_patch(Box)
Box = tr.Bbox([[0, 0], [5, 5]])
particles = random_particles(1000, ax)


def animate(i):
    dt = 0.01
    r0=0.3
    v0=3
    move_particles(particles, dt, v0, r0, Box)
    return list(map(lambda x: x.object, particles))


# create animation using the animate() function
myAnimation = animation.FuncAnimation(fig, animate, \
                                      frames=1000, blit=True, repeat=True)

myAnimation.save('vicsek.html', fps=30, extra_args=['-vcodec', 'libx264'])