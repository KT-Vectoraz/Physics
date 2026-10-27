# fluidsim

Roughly solves
```
du/dt + (u.grad)u = -grad(p)/rho + nu * lap(u) + f
div(u) = 0
dT/dt + (u.grad)T = alpha * lap(T)
f_y = g * beta * (T - T_amb)     (y is up; hot air rises)
```
## Build

```bash
gcc -std=c11 -O3 -Wall -Wextra -o fluid3d fluidsim.c util.c advection.c diffusion.c projection.c diagnostics.c setup.c -lm
```

Override grid resolution with `-DN=64`; default is `N=48`.

## Run

```bash
./fluid3d [steps]     
./fluid3d --test      
```
