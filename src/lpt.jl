
abstract type AbstractLPT end
abstract type AbstractFirstOrderLPT <: AbstractLPT end
struct FirstOrderLPT <: AbstractFirstOrderLPT end

function wavevectors3D(T::Type{<:Real}, dims, box_size=(2π, 2π, 2π))
    sample_rate = T.(2π .* dims ./ box_size)
    kx = fftfreq(dims[1], sample_rate[1])
    ky = fftfreq(dims[2], sample_rate[2])
    kz = fftfreq(dims[3], sample_rate[3])
    return kx, ky, kz
end


# this is really a reference implementation, should use rfft
function lpt(::FirstOrderLPT, δᵢₙᵢ::AbstractArray{T,3}, kv::NTuple) where T
    kxs, kys, kzs = kv
    ℱϕ⁽¹⁾ = fft(δᵢₙᵢ)
    bx, by, bz = similar(ℱϕ⁽¹⁾), similar(ℱϕ⁽¹⁾), similar(ℱϕ⁽¹⁾)
    Threads.@threads for ix in eachindex(kxs)
        kx = kxs[ix]
        kx² = kx^2
        for iy in eachindex(kys)
            ky = kys[iy]
            kx²_plus_ky² = kx² + ky^2
            for iz in eachindex(kzs)
                kz = kzs[iz]
                k² = kx²_plus_ky² + kz^2
                if k² > 0
                    im_ℱϕ⁽¹⁾_over_k² = im * ℱϕ⁽¹⁾[ix,iy,iz] / k²
                    bx[ix,iy,iz] = im_ℱϕ⁽¹⁾_over_k² * kx
                    by[ix,iy,iz] = im_ℱϕ⁽¹⁾_over_k² * ky
                    bz[ix,iy,iz] = im_ℱϕ⁽¹⁾_over_k² * kz
                else
                    bx[ix,iy,iz] = 0
                    by[ix,iy,iz] = 0
                    bz[ix,iy,iz] = 0
                end
            end
        end
    end
    return ifft(bx), ifft(by), ifft(bz)
end
