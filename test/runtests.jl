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

            @test isapprox(sdynss[:,end], [-0.535249, -0.133732, 0.811530], atol=1e-5)
        end
    
        @testset "Classical single spin steady state" begin
            sz_clgibbs(S0, T) = -(coth(S0/T) - T/S0)

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

            for S0 in [1/2, 1, 3/2, 2]
                for n in eachindex(T)
                    noise = ClassicalNoise(T[n]/S0);
                    s = zeros(navg, 3)
                    for i in 1:navg
                        bfields = [bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise),
                                bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise),
                                bfield(round(Int, 2*tend/Δt)+10, Δt/2, J, noise)];
                        sol = diffeqsolver(s0, tspan, J, bfields, matrix; S0=S0, saveat=saveat, alg=Vern7(), atol=1e-6, rtol=1e-6);
                        s[i, :] = mean(Array(sol), dims=2)
                    end
                    Sss = mean(s, dims=1)

                    if !isapprox(Sss[3], sz_clgibbs(S0, T[n]/S0), atol=1e-2)
                        println("filed with: S0 = $S0, T = $(T[n])")
                    end
                    @test isapprox(Sss[1], 0, atol=1e-2)
                    @test isapprox(Sss[2], 0, atol=1e-2)
                    @test isapprox(Sss[3], sz_clgibbs(S0, T[n]/S0), atol=1e-2)
                end
            end
        end

        @testset "Reproduce ASH paper" begin
            prm1 = LorentzianSD(10.0, 7.0, 5.0)
            prm2 = LorentzianSD(0.16, 1.4, 0.5)

            Ta = 0.07433895689641692
            Tb = 14.8677913792833860
            noise_a = QuantumNoise(Ta)
            noise_b = QuantumNoise(Tb)

            S0a = 1/2
            S0b = 200/2
        
            M = IsoCoupling(1)
            s0 = [1., 0., 0.]

            tmax = 2π*48
            Δt = 0.15
            N = floor(Int, tmax/Δt)
            tspan = (0., tmax)
            saveat = (0:1:N)*Δt
            idxavg = 1000

            navg = 256
            atol = 1e-6
            rtol = 1e-6

            function ssavg(Jlor, noise, S0)
                sols = zeros(navg)
                for i in 1:navg
                    bfields = [bfield(N+idxavg, Δt, Jlor, noise),
                               bfield(N+idxavg, Δt, Jlor, noise),
                               bfield(N+idxavg, Δt, Jlor, noise)]
                    sol = diffeqsolver(s0, tspan, Jlor, bfields, M;
                                       S0=S0, saveat=saveat, alg=Vern7(),
                                       atol=atol, rtol=rtol)
                    sols[i] = mean(sol[3, end-idxavg:end])
                end
                solavg = mean(sols)
                return solavg
            end

            @test isapprox(ssavg(prm1, noise_a, S0a), -0.24, atol=2e-2)
            @test isapprox(ssavg(prm2, noise_a, S0a), -0.28, atol=2e-2)
            @test isapprox(ssavg(prm1, noise_b, S0b), -0.85, atol=1e-2)
            @test isapprox(ssavg(prm2, noise_b, S0b), -0.85, atol=1e-2)
        end
    end
end