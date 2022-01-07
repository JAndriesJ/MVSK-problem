module data_moments
using Test

export get_means_vector,
       get_covariance_matrix,
       get_skewness_matrix,
       get_kurtosis_matrix,
       run_tests


get_means_vector(R) = sum(R, dims=1)/(m(R)-1)
get_covariance_matrix(R) = [ sum(R[:,i].*R[:,j])  for i ∈ 1:n(R), j ∈ 1:n(R)] ./ (m(R)-1)^2
get_skewness_matrix(R)   = [ sum(R[:,i].*R[:,jk[1]].*R[:,jk[2]]) for i ∈ 1:n(R), jk ∈  [get_ij_pairs(R)...]  ]./ (m(R)-1)^3
get_kurtosis_matrix(R)   = [ sum(R[:,ij[1]].*R[:,ij[2]].*R[:,kl[1]].*R[:,kl[2]]) for ij ∈ [get_ij_pairs(R)...], kl ∈  [get_ij_pairs(R)...]  ] ./ (m(R)-1)^4

## utility
m(R) = size(R)[1] ; n(R) = size(R)[2] # R ∈ ℝᵐˣⁿ
get_ij_pairs(R) = [ (j,i)  for i ∈ 1:n(R), j ∈ 1:n(R)] #

# Tests
function run_tests()
    @testset "special case" begin
        M = 100
        N = 10
        R = ones(M,N)
        means_vector      = get_means_vector(R)
        (M/(M-1))*ones(N,1)
        @test means_vector == (M/(M-1))*ones(N,1)'
        covariance_matrix = get_covariance_matrix(R)
        @test covariance_matrix == (M/(M-1)^2)*ones(N,N)
        skewness_matrix = get_skewness_matrix(R)
        @test skewness_matrix == (M/(M-1)^3)*ones(N,N^2)
        kurtosis_matrix   = get_kurtosis_matrix(R)
        @test kurtosis_matrix == (M/(M-1)^4)*ones(N^2,N^2)
    end

    @testset "moment matrices sizes" begin
        M = 100
        for N ∈ rand(7:12,6)
            R = rand(M,N)
            means_vector      = get_means_vector(R)
            @test length(means_vector) == N
            covariance_matrix = get_covariance_matrix(R)
            @test size(covariance_matrix) == (N,N)
            skewness_matrix = get_skewness_matrix(R)
            @test size(skewness_matrix) == (N,N^2)
            kurtosis_matrix   = get_kurtosis_matrix(R)
            @test size(kurtosis_matrix) == (N^2,N^2)
        end
    end

    @testset "Spec instance" begin
        R = [0.00936195  -0.0226331     0.00645874   -0.00126296    0.0133514     0.00249612 
            -0.0172047   -0.00636151    0.0157301     0.00392595   -0.0134888     0.00442713 
            -0.00233563   0.000999394   0.000228636   0.0202542     0.00634332    0.00152948 
            0.00793173  -0.0118861     0.0081441    -0.00706441   -0.00490681   -0.00392801 
            -0.053843    -0.00379494   -0.026317     -0.0127871    -0.0184033    -0.00557188 
            0.00542937  -0.00673573   -0.00508969   -0.0149126     0.00742315   -0.00107126 
            -0.0214208   -0.00800546   -0.0130998    -0.00468673   -0.0166464     0.00870028 
            -0.00391231  -0.0636361     0.000266211  -0.0139817    -0.0184761     0.00250033 
            -0.0216606   -0.0316247     0.0710702    -0.0178674    -0.00929719    0.000559572
            0.00221541  -0.0125524     0.0063527     0.00529857    0.0939015    -0.00332483 
            0.00683914  -0.00212513    0.00749073   -0.00605857    0.0350173     0.0060822  
            -0.009043     0.0153877    -0.0105436    -0.0109066    -0.0200857    -0.00688743 
            0.00415974  -0.0114596     0.00443583   -0.00295516   -0.00238862   -0.0167532  
            0.0075719   -0.00511722   -0.00304662    0.0379298    -0.0193662    -0.00445668 
            0.0131503   -0.00663526    0.0133312     0.00292408    0.0207832     0.00922465 
            0.00183876  -0.000303226  -0.00302128    0.000882107   0.00684295   -0.00677077 
            0.0137556    0.00573116    0.00963789    0.00213426    0.0325437     0.0219513  
            -0.0115745   -0.0107992    -0.00906186   -0.000765959   0.0126879    -0.00306371 
            -0.00343246  -0.0245152     0.00542973   -0.00978953    0.00482674   -0.00671045
            -0.0041411    0.000934187  -0.00483661   -0.0270348    -0.0129201    -0.0134256
            0.00729811   0.0360607    -0.0104441     0.0258354    -0.000720657   0.0123256
            -0.0149842   -0.00521257   -0.0120523    -0.00101282   -0.0234189    -0.00681499
            -9.79551e-5   0.0103868     0.00162744   -0.00580225    0.00137245   -0.00183635
            -0.00959408  -0.00486534   -0.00806891   -0.0107534    -0.000228842  -0.00691363
            0.0127866    0.00976981    0.00369384    0.00423154    0.00436418    0.0083415
            -0.0084137    0.0116491    -0.0126887     0.000895804  -0.0103714    -0.000838645
            0.0218422    0.0141859     0.0139886     0.025066      0.0089295     0.0137369
            0.00652916   0.0165108     0.000685237  -0.0699402     0.00413015    0.0095698
            -0.00225003   0.00687188    0.0124733     0.0347433     0.00749794   -0.000796934
            0.0111292   -0.00165591   -0.00168241    0.0316702    -0.00549773    0.00486077]
        means_vector      = get_means_vector(R) 
        covariance_matrix = get_covariance_matrix(R)
        skewness_matrix   = get_skewness_matrix(R)
        kurtosis_matrix   = get_kurtosis_matrix(R)
    end
end


end


