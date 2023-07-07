using SpiDy
using LinearAlgebra
using DifferentialEquations
using Test

@testset "SpiDy.jl" begin
    @testset "Spin dynamics runs" begin
        Δt = 0.1
        N = 1000
        tspan = (0, N*Δt)
        saveat = (0:1:N)*Δt
        J = LorentzianSD(10, 7, 5)
        Cw = AnisoCoupling([-sin(π/4) 0 0; 0 0 0; cos(π/4) 0 0]);
        T = 2.0
        cl_noise = ClassicalNoise(T);
        qu_noise = QuantumNoise(T);
        @test spectrum(cl_noise)(10.0) ≈ 0.4
        @test spectrum(qu_noise)(10.0) ≈ 1.0135673098126083

        nspin = 4
        s0 = zeros(3*nspin)
        for i in 1:nspin
            s0[1+(i-1)*3:3+(i-1)*3] = normalize([0.1*rand(), 0.1*rand(), -1])
        end
        J0 = 0.5
        JH = Nchain(nspin, J0)

        # only shared noise; classical
        bfields = [bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, Cw; JH=JH, saveat=saveat);
        @test sol.retcode == SciMLBase.ReturnCode.Success

        # only shared noise; quantum
        bfields = [bfield(N, Δt, J, qu_noise), bfield(N, Δt, J, qu_noise), bfield(N, Δt, J, qu_noise)];
        sol = diffeqsolver(s0, tspan, J, bfields, Cw; JH=JH, saveat=saveat);
        @test sol.retcode == SciMLBase.ReturnCode.Success

        # only local noise
        sol = diffeqsolver(s0, tspan, J, [bfields, bfields, bfields, bfields], Cw; JH=JH, saveat=saveat);
        @test sol.retcode == SciMLBase.ReturnCode.Success

        # both local and shared noise
        sol = diffeqsolver(s0, tspan, J, J, [bfields, bfields, bfields, bfields], bfields, Cw; JH=JH, saveat=saveat);
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end
end