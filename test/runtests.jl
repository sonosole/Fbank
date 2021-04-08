using Fbank
using Test

@testset "test extracting fbank features" begin
function testfeat(N::Int)
    p = fbankparams()
    getfeat = initfbank(p)
    wav = zeros(1,16000*N)

    tic = time()
    feat = getfeat(wav)
    toc = time()
    sec = (toc-tic)
    rt = floor(N / sec)
    println("===============================================")
    println("time to extrac features from 1 sec wav: ",sec," s")
    println(" rts to extrac features from 1 sec wav: ", rt)
    println("===============================================")
    return 1
end

@test testfeat(1)

end
