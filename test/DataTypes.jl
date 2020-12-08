@testset "AtomicProperties" begin
    t1 = Network_qse.AtomicProperties(2, 3, 1.0 , 10.0, x-> x^2)
    n = Network_qse.AtomicProperties(0, 1, 0.5, 8.071, x-> x^2)
    T = [1.0,2.0,3.0,4.0]
    g = [1.0,4.0,9.0,16.0]
    h  = Network_qse.AtomicProperties(2, 3, 1.0, 10.0, t-> Network_qse.LinearInterpolation(T, g)(t))
    @test t1.name == "He3"
    @test n.M*Network_qse.meverg/Network_qse.c^2 / 1.674920e-24 < 1.01
    @test t1.Z == 2
    @test t1.A == 3
    @test t1.s == 1.0
    @test t1.ω(2) == 4
    @test h.ω(2.0) == 4.0
end

@testset "ThermoProperties" begin
    @test 1 == 1
end
@testset "Func" begin
    @test 1 == 1
end

@testset "StepParameter" begin
    @test 1 == 1
end
