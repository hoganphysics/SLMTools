using Test
using SLMTools.ImageProcessing

@testset "ImageProcessing tests" begin
    # assumes the image folder testdata/test_images/ exists and has exactly three bitmaps (1216,1936) in it, named 01.bmp, 02.bmp, 03.bmp
    imgs, img_names = getImagesAndFilenames("test/test_data/test_images_A/", ".bmp")
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

    @testset "getCamGrid, sweepCoords, findCOM, and lineFit tests" begin
        # assumes the image folder testdata/test_images_B/ exists and has the images vandy put there for testing
        # getCamGrid includes a test for sweepCoords, so we don't need to test that separately
        # similarly, findCom and lineFit are tested in sweepCoords, so we don't need to test them separately
        divImgs, divImgNames = getImagesAndFilenames("./test/test_data/test_images_B/OffsetParabolas/", ".bmp")
        linImgs, linImgNames = getImagesAndFilenames("./test/test_data/test_images_B/LinearPhases/", ".bmp")

        x0, xhspan = 1050, 150 # x0 is the center of the ROI, xhspan is half the width of the ROI, use 1050, 150 for the test images
        y0, yhspan = 690, 150 # y0 is the center of the ROI, yhspan is half the height of the ROI, use 690, 150 for the test images
        roi = ((y0-yhspan):(y0+yhspan), (x0-xhspan):(x0+xhspan))
        indstart = 19
        camgrid, angle = getCamGrid(linImgs[indstart:end], [parse(Int, n[1:2]) for n in linImgNames[indstart:end]], dxcam; roi=roi)
        @test angle ≈ 0.02234140096699088    
        @test isapprox(camgrid[1], 332.8237653385899:5.86:338.6837653385899)
    end
end