from shooting_class import shoot
import math


def main():
    distance = 1.5e3
    v0 = 150.
    theta = math.pi/4.
    sh = shoot(distance,v0)
    tlist,ylist = shoot.analytical_height(theta,2.)
    print(ylist)

main()
