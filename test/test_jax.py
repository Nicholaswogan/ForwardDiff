import jax
from jax import numpy as np
from scipy.io import FortranFile

def func_operators(x):
    res = x + x + 3.0 + x
    res = res - x - 3.0 - x
    res = res*x + 2.0*x + x*-5.0
    res = res/x + 2.0/x + x/5.0
    res = res**x + res**1.5
    return res

def func_intrinsics1(x):
    res = np.abs(x)
    res = res + np.cos(x)
    res = res + np.exp(x)
    res = res + np.log(x)
    res = res + np.log10(x)
    res = res + np.sin(x)
    res = res + np.tan(x)
    res = res + np.sqrt(x)
    return res

def func_intrinsics2(x):
    res = np.arccos(x)
    res = res + np.arcsin(x)
    res = res + np.arctan(x)
    res = np.maximum(res, x)
    res = np.maximum(res, 1.0)
    res = np.maximum(1.0, res)
    res = np.minimum(res, res)
    res = np.minimum(res, 2.0)
    res = np.minimum(2.0, res)
    return res

def func_grad1(x):
    res = x[0]*x[0]*x[1] + x[0] + x[1]
    return res

def func_grad2(x):
    res = np.sum(x*3.14)
    return res

def test():
    fil = FortranFile('test.dat','r')

    x = np.array(3.0,dtype=np.float32)
    f = func_operators(x)
    dfdx = jax.grad(func_operators)(x)
    f1, dfdx1 = fil.read_record(np.float64)
    print(f/f1,dfdx/dfdx1)
    assert np.isclose(f,f1) and np.isclose(dfdx,dfdx1)

    x = np.array(2.0,dtype=np.float32)
    f = func_intrinsics1(x)
    dfdx = jax.grad(func_intrinsics1)(x)
    f1, dfdx1 = fil.read_record(np.float64)
    print(f/f1,dfdx/dfdx1)
    assert np.isclose(f,f1) and np.isclose(dfdx,dfdx1)

    x = np.array(0.1,dtype=np.float32)
    f = func_intrinsics2(x)
    dfdx = jax.grad(func_intrinsics2)(x)
    f1, dfdx1 = fil.read_record(np.float64)
    print(f/f1,dfdx/dfdx1)
    assert np.isclose(f,f1) and np.isclose(dfdx,dfdx1)

    x = np.array([1, 2],dtype=np.float32)
    f = func_grad1(x)
    dfdx = jax.grad(func_grad1)(x)
    tmp = fil.read_record(np.float64)
    f1, dfdx1 = tmp[0],tmp[1:]
    print(f/f1,dfdx/dfdx1)
    assert np.isclose(f,f1) and np.all(np.isclose(dfdx,dfdx1))

    x = np.array([3, 4],dtype=np.float32)
    f = func_grad2(x)
    dfdx = jax.grad(func_grad2)(x)
    tmp = fil.read_record(np.float64)
    f1, dfdx1 = tmp[0],tmp[1:]
    print(f/f1,dfdx/dfdx1)
    assert np.isclose(f,f1) and np.all(np.isclose(dfdx,dfdx1))

    fil.close()

if __name__ == "__main__":
    test()
