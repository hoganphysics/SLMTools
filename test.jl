using SLMTools
using Test

saveAs8BitBMP(convert.(Int64, floor.(rand(100,100).*255)), "test.bmp")