using SLMTools, FileIO
# cd("z:\\home\\Vandy\\code\\julia\\SLMTools\\examples")
divImgs, divImgNames = getImagesAndFilenames("./test/test_data/test_images_B/OffsetParabolas/", ".bmp");
linImgs, linImgNames = getImagesAndFilenames("./test/test_data/test_images_B/LinearPhases/", ".bmp");

x0, xhspan = 1050, 150
y0, yhspan = 690, 150
roi = ((y0-yhspan):(y0+yhspan), (x0-xhspan):(x0+xhspan))

# indstart = 19
# camgrid, angle = getCamGrid(linImgs[indstart:end], [parse(Int, n[1:2]) for n in linImgNames[indstart:end]], dxcam; roi=roi)
angle = 0.02234140096699088    
camgrid= (332.8237653385899:5.86:338.6837653385899,)
interps = dualate(divImgs[end], Lslm, camgrid, θ, flambda; roi=roi)