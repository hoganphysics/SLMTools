using Test
using SLMTools.ImageProcessing

@testset "ImageProcessing tests" begin
    # assumes the image folder testdata/test_images/ exists and has exactly three bitmaps (1216,1936) in it, named 01.bmp, 02.bmp, 03.bmp
    imgs, img_names = getImagesAndFilenames(raw"test_data\test_images_A\\", ".bmp")
    
    @testset "getImagesAndFilenames tests" begin
        @test size(imgs) == (3,)
        @test all(size(img) == (1216, 1936) for img in imgs[1:3])
        @test img_names == ["01.bmp", "02.bmp", "03.bmp"]
    end

    fa1 = imageToFloatArray(imgs[1])   

    @testset "imageToFloatArray tests" begin
        @test size(fa1) == (1216, 1936)
        @test typeof(fa1) == Matrix{Float64}
        @test itfa == imageToFloatArray
    end
end