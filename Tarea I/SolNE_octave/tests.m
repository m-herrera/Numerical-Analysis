pkg load symbolic
[x1, k1] = sne_ud_1(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)
[x2, k2] = sne_ud_2(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x3, k3] = sne_ud_3(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x4, k4] = sne_ud_4(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x5, k5] = sne_ud_5(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x6, k6] = sne_ud_6(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x7, k7] = sne_ud_7(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)

[x1, k1] = sne_fd_1(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x2, k2] = sne_fd_2(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x3, k3] = sne_fd_3(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x4, k4] = sne_fd_4(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x5, k5] = sne_fd_5(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x6, k6] = sne_fd_6(100, 0.0001, "@(x)x**2+3*x-10", 1)
[x7, k7] = sne_fd_7(100, 1, 0.0001, "@(x)x**2+3*x-10", 1)