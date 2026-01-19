import numpy as np
from numpy.polynomial import Polynomial
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline


def poly_on_grid(x, x0=0, coeffs=[1,2,3]):
    # polynomials and derivatives
    p = Polynomial(coeffs)
    dp = p.deriv()

    # polynmial and derivative on grid
    return p(x-x0), dp(x-x0)


# product potental energy surface, 2d
def product_poly_2d(x,y, x0=0,y0=0, coeffs_x=[1,2,3], coeffs_y=[1,2,3]):

    # polynomials and derivatives
    fx, dfx = poly_on_grid(x, x0=x0, coeffs=coeffs_x)
    fy, dfy = poly_on_grid(y, x0=y0, coeffs=coeffs_y)

    return {'f': fx[:, None] @ fy[None, :], 
            'df_dx': dfx[:, None] @ fy[None, :],
            'df_dy': fx[:, None] @ dfy[None, :]}

if __name__ == "__main__":
    x= np.linspace(-1,1, 5)
    y= np.linspace(-1,1, 5)

    ret = product_poly_2d(x,y, x0=0,y0=0, 
                          coeffs_x=[1,2,3,4], 
                          coeffs_y=[1,2,3])
    
    # spline interpolation on 2d grid
    Z = ret.get('f')
    spl = RectBivariateSpline(x, y, Z, kx=3, ky=3)

    X, Y = np.meshgrid(x, y, indexing="ij")

    # --- Plot ---
    fig = plt.figure(figsize=(12, 4))

    # 3D surface
    
    Z = ret.get('f')
    ax1  = fig.add_subplot(1, 3, 1, projection="3d")
    ax1.plot_surface(X, Y, Z, cmap="viridis")
    ax1.set_title("Z = f(x) · g(y)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")

    Z = ret.get('df_dx')
    ax2  = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.plot_surface(X, Y, Z, cmap="viridis")
    ax2.set_title("Z = df(x)/dx · g(y)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")

    Z = ret.get('df_dy')
    ax2  = fig.add_subplot(1, 3, 3, projection="3d")
    ax2.plot_surface(X, Y, Z, cmap="viridis")
    ax2.set_title("Z = f(x) · dg(y)/dy")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")


    # plot the thing on finer grid
    xf= np.linspace(x.min(),x.max(), 200)
    yf= np.linspace(y.min(),y.max(), 200)


    retf = product_poly_2d(xf,yf, x0=0,y0=0, 
                          coeffs_x=[1,2,3,4], 
                          coeffs_y=[1,2,3])


    Xf, Yf = np.meshgrid(xf, yf, indexing="ij")

    fig = plt.figure(figsize=(12, 4))

    # 3D surface
    Zf   = spl(xf, yf) 
    ax1  = fig.add_subplot(1, 3, 1, projection="3d")
    ax1.plot_surface(Xf, Yf, Zf, cmap="viridis")
    ax1.set_title("Z = f(x) · g(y)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")

    Zf = spl(xf, yf, dx=1, dy=0)    
    ax2  = fig.add_subplot(1, 3, 2, projection="3d")
    ax2.plot_surface(Xf, Yf, Zf, cmap="viridis")
    ax2.set_title("Z = df(x)/dx · g(y)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")

    Zf = spl(xf, yf, dx=0, dy=1) 
    ax2  = fig.add_subplot(1, 3, 3, projection="3d")
    ax2.plot_surface(Xf, Yf, Zf, cmap="viridis")
    ax2.set_title("Z = f(x) · dg(y)/dy")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")

    plt.show()