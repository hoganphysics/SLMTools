using SLMTools, FileIO
using Images: Gray, RGB
divImgs, divImgNames = getImagesAndFilenames("./test/test_data/test_images_B/OffsetParabolas/", ".bmp");
linImgs, linImgNames = getImagesAndFilenames("./test/test_data/test_images_B/LinearPhases/", ".bmp");

x0, xhspan = 1050, 150
y0, yhspan = 690, 150
roi = ((y0-yhspan):(y0+yhspan), (x0-xhspan):(x0+xhspan))
# roi = (550:950, 825:1225)
# [d[roi...] for d in divImgs] # display the cropped image, we want the zero order spot to be cropped out, with the first order spot visible
indstart = 19
camgrid, angle = getCamGrid(linImgs[indstart:end], [parse(Int, n[1:2]) for n in linImgNames[indstart:end]], dxcam; roi=roi)

# angle = 0.02234140096699088
# camgrid = (332.8237653385899:5.86:338.6837653385899,)
# θ= angle
# camGrid = camgrid 
# imgs = divImgs[end-2:end]
# if !isnothing(roi)
#     if size(imgs[1]) != length.(roi)
#         imgs = [i[roi...] for i in imgs]
#     end
#     if length.(camGrid) != length.(roi)
#         camGrid = ((camGrid[i][roi[i]] for i = 1:length(camGrid))...,)
#     end
# end
# size(imgs[1]) == length.(camGrid) || error("Camera grid and image size unequal.")
# interps = [linear_interpolation(camGrid, imgs[j], extrapolation_bc=zero(typeof(imgs[1][1]))) for j in 1:length(imgs)]
# dL = dualShiftLattice(slmGrid, flambda)
# Lx0, Lxs, Ly0, Lys = dL[1][1], step(dL[1]), dL[2][1], step(dL[2])
# R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
# r0, dx, dy = (R * [Lx0, Ly0]), Lxs * R[:, 1], Lys * R[:, 2]

# [[interp((r0 + dx * I[1] + dy * I[2])...) for I in CartesianIndices(length.(dL))] for interp in interps]