using SpiDy
using LinearAlgebra
using DifferentialEquations
using Statistics
using PolynomialRoots
using Random
using Test

@testset "SpiDy.jl" begin
    @testset "Local/shared bath syntax tests" begin
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
        Jlist = [J, J, J, J, J]
        bfieldlist = [bfields, bfields, bfields, bfields, bfields]
        bcoupling = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,1,1,1]]
        matrix = [Cw, Cw, Cw, Cw, Cw]
        sol = diffeqsolver(s0, tspan, Jlist, bfieldlist, bcoupling, matrix; JH=JH, saveat=saveat);
        @test sol.retcode == SciMLBase.ReturnCode.Success
    end

    @testset "Spin dynamics tests" begin
        @testset "Classical T=0 test vs analytics" begin
            function compare_T0_ss(α, θ0)
                Δt = 0.1
                N = 80_000
                tspan = (0, N*Δt)
                saveat = (9*N÷10:1:N)*Δt
            
                J = LorentzianSD(α, 7, 5)
                Cw = AnisoCoupling([-sin(θ0) 0 0; 0 0 0; cos(θ0) 0 0]);
                T = 0.0
                cl_noise = ClassicalNoise(T);
                bfields = [bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise)];
            
                s0 = normalize(rand(3))
                sol = diffeqsolver(s0, tspan, J, bfields, Cw; saveat=saveat, alg=Tsit5(), atol=1e-5, rtol=1e-5);
                sdynss = mean(Array(sol), dims=2)
            
                # analytic solution of T=0 cl steady-state
                norm_eq_coeffs(ζ, θ0) = [-1, 4*ζ, 1 - 4*ζ^2, -4*ζ*sin(θ0)^2, 4*ζ^2*sin(θ0)^2]
                sx_ss(Z, ζ, θ0) = 2*Z^2*ζ*sin(θ0)*cos(θ0)/(1 - 2*Z*ζ)
                sz_ss(Z, ζ, θ0) = Z + 2*Z^2*ζ*cos(θ0)^2/(1 - 2*Z*ζ)
                function s_ss(ζ, θ0)
                    Z = real(roots(norm_eq_coeffs(ζ, θ0))[3])
                    [sx_ss(Z, ζ, θ0), 0,  sz_ss(Z, ζ, θ0)]
                end
                sanalytss = s_ss(SpiDy.reorganisation_energy(J), π/2 - θ0)
            
                return sdynss, sanalytss
            end 

            for (α,θ0) in [(10, π/4), (20, 3*π/8), (17, 5*π/9)]
                sdynss, sanalytss = compare_T0_ss(α, θ0)
                @test norm(sdynss - sanalytss) < 1e-3
            end
        end

        @testset "Stochastic test with fixed seed" begin
            Δt = 0.1
            N = 100
            tspan = (0, N*Δt)
            saveat = (9*N÷10:1:N)*Δt
        
            J = LorentzianSD(10, 7, 5)
            Cw = AnisoCoupling([-sin(π/4) 0 0; 0 0 0; cos(π/4) 0 0]);
            T = 0.1
            cl_noise = QuantumNoise(T);

            Random.seed!(1398465)
            bfields = [bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise), bfield(N, Δt, J, cl_noise)];
        
            s0 = normalize(rand(3))
            sol = diffeqsolver(s0, tspan, J, bfields, Cw; S0=1, saveat=saveat, alg=Vern7(), atol=1e-8, rtol=1e-8);
            sdynss = mean(Array(sol), dims=2)

            @test isapprox(sdynss[:,end], [-0.535859, -0.109671, 0.814815], atol=1e-5)
        end
    
        @testset "Classical single spin steady state" begin
            S0 = 1/2
            sz_clgibbs(T) = -(coth(S0/T) - T/S0)

            tstr, tend = 50, 200
            Δt = 0.5
            tspan = (0., tend)
            saveat = tstr:Δt:tend

            α, ω0, Γ = 1, 2, 3
            J = LorentzianSD(α, ω0, Γ);
            matrix = IsoCoupling(1);

            T = 10 .^ LinRange(-2, 2, 5);

            navg = 200
            s0 = normalize([0.1, 0.0, -1.0]);

            # although the results being tested here should be independent of
            # the random seed, we set it to a fixed value for reproducibility,
            # especially given the low number of samples used to speed up the
            # testing
            Random.seed!(2827465)

            for n in eachindex(T)
                noise = ClassicalNoise(T[n]);
                s = zeros(navg, 3)
                for i in 1:navg
                    bfields = [bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise),
                               bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise),
                               bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise)];
                    sol = diffeqsolver(s0, tspan, J, bfields, matrix; saveat=saveat, alg=Vern7(), atol=1e-6, rtol=1e-6);
                    s[i, :] = mean(Array(sol), dims=2)
                end
                Sss = mean(s, dims=1)

                @test isapprox(Sss[1], 0, atol=1e-2)
                @test isapprox(Sss[2], 0, atol=1e-2)
                @test isapprox(Sss[3], sz_clgibbs(T[n]), atol=1e-2)
            end
        end
    end
end