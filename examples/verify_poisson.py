"""Week02 独立 SI Poisson MMS：直线/二次/正弦，及 SI↔nd 行尺度。

此教学脚本使用 SI [C/m³,V] 混合残差，绝不作为 W3 生产 F(X)。
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

from minitcad.core.constants import q, k_b_j_per_k, epsilon_0_f_per_m
from minitcad.discretization.poisson import poisson_residual_nd
from examples.poisson_si_oracle import poisson_residual_si_mms


def verify() -> tuple[float,...]:
    epsilon = 11.7 * epsilon_0_f_per_m
    L = 1e-6
    x = np.array([0.0,0.07,0.24,0.63,1.0])*L
    linear = 0.2 + 0.7*x/L
    err = poisson_residual_si_mms(x,linear,np.zeros_like(x),epsilon,linear[0],linear[-1])
    assert np.max(abs(err)) < 1e-7
    quadratic = -0.3 + 0.4*x/L + 0.8*(x/L)**2
    rho = np.full(x.size,-1.6*epsilon/L**2)
    err = poisson_residual_si_mms(x,quadratic,rho,epsilon,quadratic[0],quadratic[-1])
    assert np.max(abs(err)) < 1e-7

    errors = []
    for m in (17,33,65,129):
        grid = np.linspace(0.0,L,m)
        dx = L/(m-1)
        exact = np.sin(np.pi * grid/L)
        density = epsilon*(np.pi/L)**2*exact
        dim = m-2
        a = diags([
            np.full(dim-1,-epsilon/dx**2),
            np.full(dim,2*epsilon/dx**2),
            np.full(dim-1,-epsilon/dx**2)
        ],[-1,0,1],format='csc')
        calculated = np.zeros(m)
        calculated[1:-1] = spsolve(a,density[1:-1])
        errors.append(float(np.max(abs(calculated-exact))))
        residual = poisson_residual_si_mms(grid,calculated,density,epsilon,0.0,0.0)
        assert np.max(abs(residual[1:-1]))/max(1,np.max(abs(density))) < 1e-11
    ratios = [errors[i]/errors[i+1] for i in range(len(errors)-1)]
    assert all(3.8 <= ratio <= 4.2 for ratio in ratios)

    T,ni_cm3 = 300.0,1e10
    VT = k_b_j_per_k*T/q
    l0_cm = np.sqrt((epsilon/100)*VT/(q*ni_cm3))
    xi = np.array([0.0,0.25,0.75,1.5])
    u = 0.5+0.2*xi-0.3*xi**2
    n,p = np.array([1.1,1.2,2.3,2.4]),np.array([1.4,1.2,0.8,0.6])
    cbar = np.array([-0.4,0.7])
    fnd = poisson_residual_nd(xi,u,n,p,cbar,u[0]+0.12,u[-1]-0.13)
    charge = q*(ni_cm3*1e6)
    rho = charge*(np.array([0,*cbar,0])+p-n)
    fsi = poisson_residual_si_mms(
        xi*l0_cm/100,u*VT,rho,epsilon,(u[0]+0.12)*VT,(u[-1]-0.13)*VT)
    np.testing.assert_allclose(fsi[1:-1]/(-charge),fnd[1:-1],atol=1e-12,rtol=1e-12)
    np.testing.assert_allclose(fsi[[0,-1]]/VT,fnd[[0,-1]],atol=1e-12,rtol=1e-12)
    return tuple(ratios)


if __name__ == '__main__':
    ratios = verify()
    print('SI MMS & ND scaling: PASS')
    print('Sine refinement error ratios (17/33/65/129 nodes):',
          ', '.join(f'{value:.8f}' for value in ratios))
