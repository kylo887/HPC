#include <stdio.h>
#include <stdlib.h>
int main()
{
    unsigned seed;
    int r,a,b;
    printf("input:");
    scanf("%u", &seed);
    srand(seed);
    for(a=0;a<20;a++)
    {
        for(b=0;b<5;b++)
        {
            r=rand();
            printf("%dt",r);
        }
        putchar('n');
    }
    return(0); 
}

/*
n = Number of particles in simulation;
for i = 0...n;

ax = 0.0;
ay = 0.0;
az = 0.0;
for j=0...n do
∆x = x[j]-x[i];
∆y = y[j]-y[i];
∆z = z[j]-z[i];
invr = 1:0=p(∆x2 + ∆y2 + ∆z2 + eps);
invr3 = invr3;
f=m[j]*invr3;
ax += f*∆x;
ay += f*∆y;
az += f*∆x;
end
xnew[i] = x[i] + dt*vx[i] + 0.5*dt*dt*ax;
ynew[i] = y[i] + dt*vy[i] + 0.5*dt*dt*ay;
znew[i] = z[i] + dt*vz[i] + 0.5*dt*dt*az;
vx[i] += dt*ax;
vy[i] += dt*ay;
vz[i] += dt*az;
end
for i = 0...n do
x[i] = xnew[i];
y[i] = ynew[i];
z[i] = znew[i];
end
*/
