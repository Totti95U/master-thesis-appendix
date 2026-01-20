
function matrix2str(A::Matrix, align::Symbol=:left)
    w_each_col = maximum(length.(string.(A)), dims=1)
    str = ""
    for i in axes(A, 1)
        if i > 1
            str *= (i == size(A, 1)) ? "⎣" : "⎢"
        else
            str *= (size(A, 1) == 1) ? "[" : "⎡"
        end

        for j in axes(A, 2)
            if align == :left
                str *= string(A[i, j])
                str *= " "^(w_each_col[j] - length(string(A[i, j])))
            else
                str *= " "^(w_each_col[j] - length(string(A[i, j])))
                str *= string(A[i, j])
            end

            if j < size(A, 2)
                str *= " "
            end
        end

        if i < size(A, 1)
            str *= (i == 1) ? "⎤" : "⎥"
            str *= "\n"
        else
            str *= (size(A, 1) == 1) ? "]" : "⎦"
        end
    end
    return str
end

function _pos2brailleΔcode(i, j)
    if i <= 3
        return 2^(i-1 + 3*(j-1))
    else
        return 2^(6 + (j-1))
    end
end

"""
    matrix_to_braille(matrix::Matrix{T}, [maxHeight], [maxWidth]) where T <: Number

Convert a matrix to a string representation using braille characters.
Each braille character represents a 2×4 block of the matrix, where non-zero elements are shown as dots.
If the matrix is too large, it will be scaled down to fit within `maxWidth` and `maxHeight`.

# Examples
```repl
julia> println(matrix2braille(Matrix(I, 16, 16)))
⎡⠑⢄⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠑⢄⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠑⢄⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠑⢄⎦
```
"""
function matrix2braille(A::Matrix, maxHeight::Union{Int, Nothing}=nothing, maxWidth::Union{Int, Nothing}=nothing)
    # This function originally written in SparseArrays.jl
    # We wrote a slightly modified version for our purpose here

    empty_braille_code = 10240

    if !isnothing(maxWidth)
        maxWidth = max(maxWidth-2, 1)# for brackets
    end

    h, w = size(A)
    local scaleWidth::Int, scaleHeight::Int
    if isnothing(maxHeight)
        scaleHeight = h
        if !isnothing(maxWidth) && (w > 2maxWidth)
            s = 2maxWidth / w
            scaleWidth = 2maxWidth
            scaleHeight = floor(Int, s * h)
        else
            scaleWidth = w
        end
    elseif isnothing(maxWidth)
        scaleWidth = w
        if !isnothing(maxHeight) && (h > 4maxHeight)
            s = 4maxHeight / h
            scaleHeight = 4maxHeight
            scaleWidth = floor(Int, s * w)
        else
            scaleHeight = h
        end
    elseif h > 4maxHeight || w > 2maxWidth
        s = min(4maxHeight / h, 2maxWidth / w)
        scaleWidth = floor(Int, s * w)
        scaleHeight = floor(Int, s * h)
    else
        scaleWidth = w
        scaleHeight = h
    end

    scaleHeight = max(scaleHeight, 1)
    scaleWidth = max(scaleWidth, 1)

    local brailleGrid::Matrix{UInt16}
    not_row_vec = (scaleHeight - 1)÷4 + 1 > 1

    if not_row_vec
        brailleGrid = fill(UInt16(empty_braille_code), (scaleWidth-1)÷2 + 4, (scaleHeight - 1)÷4 + 1)
        brailleGrid[1, 1] = '⎡'
        brailleGrid[1, 2:end-1] .= '⎢'
        brailleGrid[1, end] = '⎣'
        brailleGrid[end-1, 1] = '⎤'
        brailleGrid[end-1, 2:end-1] .= '⎥'
        brailleGrid[end-1, end] = '⎦'
        brailleGrid[end, 1:end-1] .= '\n'
    else
        brailleGrid = fill(Char(empty_braille_code), (scaleWidth-1)÷2 + 3, 1)
        brailleGrid[1, 1] = '['
        brailleGrid[end, 1] = ']'
    end

    rowscale = max(1, scaleHeight - 1) / max(1, size(A, 1) - 1)
    colscale = max(1, scaleWidth - 1) / max(1, size(A, 2) - 1)

    for j in axes(A, 2)
        sj = round(Int, (j-1) * colscale + 1)
        for i in axes(A, 1)
            si = round(Int, (i-1) * rowscale + 1)

            k = (sj - 1) ÷ 2 + 2
            l = (si - 1) ÷ 4 + 1
            brailleGrid[k, l] |= (A[i, j] != 0) ? _pos2brailleΔcode((si - 1) % 4 + 1, (sj - 1) % 2 + 1) : 0
        end
    end

    if not_row_vec
        reshape(brailleGrid, 1, :)
        return join(Char.(brailleGrid[1:end-1]), "")
    else
        return join(Char.(brailleGrid), "")
    end
end

#=
⎡⠑⢄⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠀⠑⢄⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠑⢄⠀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠑⢄⎦
=#

function _my_sparse_rand(m, n, p=0.3)
    A = zeros(m, n)
    for i in 1:m, j in 1:n
        A[i, j] = rand() < p ? 1 : 0
    end
    return A
end
