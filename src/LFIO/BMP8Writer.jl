# Requires:
# using ColorTypes, FixedPointNumbers   # Gray, N0f8

# ── Size constants (bytes) ──────────────────────────────────────────
const FILE_HDR  = 14                # BMPFILEHEADER
const DIB_HDR   = 40                # BITMAPINFOHEADER
const PALETTE   = 256 * 4           # 256 entries × BGRA
const HEADER_SZ = FILE_HDR + DIB_HDR + PALETTE

export save_gray8bmp

"""
    save_gray8bmp(path, img)

Write `img` (2‑D, convertible to `Gray{N0f8}`) to `path`
as an uncompressed **8‑bit BMP** (standard greyscale palette).
"""
function save_gray8bmp(path::AbstractString,
                       img::AbstractMatrix{<:ColorTypes.Gray})

    g = N0f8.(img)
    H, W   = size(g)
    rowpad = UInt8(mod(-W, 4))      # 0–3 pad bytes per row
    rowsz  = W + rowpad
    datasz = rowsz * H
    filesize = HEADER_SZ + datasz

    open(path, "w") do io
        # ↓↓↓ FILE HEADER (14 bytes) ↓↓↓
        write(io, "BM")
        write(io, UInt32(filesize))
        write(io, UInt32(0))                    # reserved
        write(io, UInt32(HEADER_SZ))            # pixel offset

        # ↓↓↓ DIB HEADER (BITMAPINFOHEADER) ↓↓↓
        write(io, UInt32(DIB_HDR))              # biSize
        write(io, Int32(W))                     # width
        write(io, Int32(H))                     # height (bottom‑up)
        write(io, UInt16(1))                    # planes
        write(io, UInt16(8))                    # bitcount
        write(io, UInt32(0))                    # BI_RGB (no compression)
        write(io, UInt32(datasz))               # image size
        write(io, UInt32(0))                    # x‑pixels per metre
        write(io, UInt32(0))                    # y‑pixels per metre
        write(io, UInt32(0))                    # colours used
        write(io, UInt32(0))                    # important colours

        # ↓↓↓ 256‑entry greyscale palette ↓↓↓
        for i in 0:255
            write(io, UInt8(i))                 # Blue
            write(io, UInt8(i))                 # Green
            write(io, UInt8(i))                 # Red
            write(io, UInt8(0))                 # Reserved
        end

        # ↓↓↓ Pixel data (bottom‑up, padded) ↓↓↓
        pad = fill(UInt8(0), rowpad)
        for r = H:-1:1
            write(io, g[r, :])
            write(io, pad)
        end
    end
    return nothing
end
