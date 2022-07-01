function diffEqSolver(N, Δt, J::GenericSD, noise::Noise, distro=nothing)
    
end



# def M_sim(barTv, S0, prm, Nsam, nsx, nsy, nsz, whatNoise, matrix):
#     w0, Gamma, A = prm
#     matrix_c_w_2 = matrix @ np.transpose(matrix)
#     bx_int = intpl.interp1d(tva_extended,
#                                 b(barTv, S0, prm, Nsam, whatNoise, nsx)[:, 0],
#                                 bounds_error='False',
#                                 fill_value=0.0)
#     by_int = intpl.interp1d(tva_extended,
#                                 b(barTv, S0, prm, Nsam, whatNoise, nsy)[:, 0],
#                                 bounds_error='False',
#                                 fill_value=0.0)
#     bz_int = intpl.interp1d(tva_extended,
#                                 b(barTv, S0, prm, Nsam, whatNoise, nsz)[:, 0],
#                                 bounds_error='False',
#                                 fill_value=0.0)
#     bn = lambda t: matrix @ np.array([bx_int(t), by_int(t), bz_int(t)])
#     def system(t, V, res):
#         s = V[0:3]
#         p = V[3:6]
#         x = V[6:9]
#         Beff = np.sign(gam) * (Bn + matrix_c_w_2 @ x + bn(t))
#         res[0] = s[1] * Beff[2] - s[2] * Beff[1]
#         res[1] = s[2] * Beff[0] - s[0] * Beff[2]
#         res[2] = s[0] * Beff[1] - s[1] * Beff[0]
#         res[3:6] = -(w0**2) * x - Gamma * p + np.sign(gam) * A * s
#         res[6:9] = p
#     resa = odeint(system, tva, inspi, rtol=ode_tol, atol=ode_tol)
#     M = np.sign(gam) * resa.values.y[:, 0:3]
#     return M