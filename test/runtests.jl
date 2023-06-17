using ComminWeath
using Test

@testset "Measured Data" begin

    @test Tanaka2015() == Tanaka2015()
    @test TaylorI() == TaylorI()
    @test TaylorIII() == TaylorIII()
    @test TaylorIV() == TaylorIV()

end

@testset "Custom structs" begin
    
    graintest = Grain()
    graintest.K = 2.
    @test graintest.K == 2.
    @test graintest.rfr == 7.

    rindtest = Rind()
    rindtest.evolve=false
    @test rindtest.evolve == false
    @test rindtest.p == 0.4

    wxtest = WxAuth()
    wxtest.k = 0
    @test iszero(wxtest.k)
    @test isone(wxtest.k_power)

    detritaltest = Detrital()
    detritaltest.cU = 100
    @test detritaltest.cU == 100.
    @test isone(detritaltest.r48)

end

@testset "ComminWeath Maths" begin
    
    linetest = linreg([1,2,3],[2,5,8])
    @test linetest.m ≈ 3
    @test linetest.b ≈ -1
    @test linetest.r² ≈ 1
    
    linetest = linreg([-1,2,3],[2,5,7])
    @test linetest.m ≈ 1.192307692307692
    @test linetest.b ≈ 3.076923076923077
    @test linetest.r² ≈ 0.9726720647773279

    cwtest = comminweath(30,Grain(),Detrital(),WxAuth(),Rind(),timeseries=0:1e5:4e5)
    @test cwtest.A234 ≈ [1.0000000000000002, 1.0242593098674027, 1.041486868925901, 1.053692051266263, 1.0623101989881465]
    @test cwtest.A230 ≈ [1.0, 1.1077301896401268, 1.123763817763784, 1.130056498093801, 1.1339556060405112] 
    @test cwtest.cU ≈ [200.0, 200.701300006131, 201.40197163645126, 202.1020154539954, 202.80143202129346]

    cwtest = comminweath(30,Grain(),Detrital(),WxAuth(),Rind())
    testdate = drawdate(1000,cwtest)
    @test testdate.A234 ≈ 1.0777523364628807
    @test testdate.A230 ≈ 1.1389189406897795
    @test testdate.cU ≈ 206.981674931185

    testdates = drawdates([500,1000],cwtest)
    @test testdates.A234 ≈ [1.0638377663360574, 1.0777523364628807]
    @test testdates.A230 ≈ [1.1332054610159914, 1.1389189406897795]
    @test testdates.cU ≈ [203.49865706542445, 206.981674931185]
    
end

@testset "Visualization Prep" begin
    
    @test calcslopes([100,600,1200]) ≈ [0.2906483540712149, 0.5209640753756838, 0.595234494683003]

    @test calcU([10,40],[500,1500]) ≈ [210.48075220579932 231.0291306235513; 202.62156467709448 207.83819039321148]

end