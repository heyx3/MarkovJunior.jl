struct CellType
    code::UInt8
    name::String
    char::Char
    color::v3f # Linear color
end

"The number of possible cell types is limited, so that some operations can be optimized into bit math"
const N_CELL_TYPES = 16

const CELL_TYPES::NTuple{N_CELL_TYPES, CellType} = tuple(
    CellType(0, "Black",   'b', v3f(0, 0, 0)),
    CellType(1, "Gray",    'g', v3f(0.5, 0.5, 0.5)),
    CellType(2, "White",   'w', v3f(1, 1, 1)),
    CellType(3, "Red",     'R', v3f(1, 0, 0)),
    CellType(4, "Green",   'G', v3f(0, 1, 0)),
    CellType(5, "Blue",    'B', v3f(0, 0, 1)),
    CellType(6, "Yellow",  'Y', v3f(1, 1, 0)),
    CellType(7, "Magenta", 'M', v3f(1, 0, 1)),
    CellType(8, "Teal",    'T', v3f(0, 1, 1)),
    CellType(9, "Orange",  'O', v3f(1, 0.5, 0)),
    CellType(10, "Pink",   'P', v3f(1, 0, 0.5)),
    CellType(11, "Olive",  'L', v3f(0, 0.5, 0.2)),
    CellType(12, "Indigo", 'I', v3f(0, 0.2, 0.5)),
    CellType(13, "Brown",  'N', v3f(0.5, 0.2, 0)),
    CellType(14, "Beige",  'E', v3f(1, 0.9, 0.8)),
    CellType(15, "SkyBlue", 'S', v3f(0.7, 0.85, 1))
)


####################################
##   Set of Cell Types

@assert(N_CELL_TYPES == 16, "Have to change how we do cell sets")

cell_type_to_bitmask(t::UInt8)::UInt16 = (UInt16(1) << UInt16(t))

"
A deterministic and optimized set of the different cell types.
You can even remove elements from it *while* iterating through it!
"
struct CellTypeSet <: Base.AbstractSet{UInt8}
    bitfield::UInt16

    CellTypeSet(bits::UInt16) = new(bits)
    CellTypeSet() = new(UInt16(0))

    CellTypeSet(types::UInt8...) = CellTypeSet(types)
    CellTypeSet(types::Tuple{Vararg{UInt8}}) = new(reduce(|, map(cell_type_to_bitmask, types), init=UInt16(0)))
    CellTypeSet(types::AbstractArray{UInt8}) = new(reduce(|, map(cell_type_to_bitmask, types), init=UInt16(0)))
end

Base.in(type::UInt8, set::CellTypeSet) = (set.bitfield & cell_type_to_bitmask(type)) != 0
Base.hasfastin(::Type{CellTypeSet}) = true

Base.IteratorSize(::Type{CellTypeSet}) = Base.HasLength()
Base.length(s::CellTypeSet) = count_ones(s.bitfield)
Base.isempty(s::CellTypeSet) = iszero(s.bitfield)

function cell_set_index_of(s::CellTypeSet, value::UInt8)::Optional{Int}
    for (i, v) in enumerate(s)
        if v == value
            return i
        end
    end
    return nothing
end

Base.empty(::CellTypeSet) = CellTypeSet()
Base.empty(::Type{CellTypeSet}, ::Type{UInt8}=UInt8) = CellTypeSet()

Base.copy(s::CellTypeSet) = s
push(s::CellTypeSet, type::UInt8) = CellTypeSet(s.bitfield | cell_type_to_bitmask(type))
delete(s::CellTypeSet, type::UInt8) = CellTypeSet(s.bitfield & (~cell_type_to_bitmask(type)))

function Base.setdiff(s::CellTypeSet, elements...)
    for iter in elements
        for type in iter
            s = delete(s, convert(UInt8, type))
        end
    end
    return s
end
function Base.intersect(s::CellTypeSet, elements...)
    return reduce(elements, init=s) do s,iter
        # Convert to our accelerated set type.
        iter_set::CellTypeSet = if iter isa CellTypeSet
            iter
        else
            s2 = CellTypeSet()
            for type in iter
                s2 = push(s2, convert(UInt8, type))
            end
            s2
        end

        # Remove anything from the output that isn't in this new set.
        for type in s
            if !in(type, iter_set)
                s = delete(s, type)
            end
        end

        return s
    end
end
function Base.union(s::CellTypeSet, elements...)
    for iter in elements
        for type in iter
            s = push(s, convert(UInt8, type))
        end
    end
    return s
end
Base.union(a::CellTypeSet, b::CellTypeSet, rest...) = union(CellTypeSet(a.bitfield | b.bitfield), rest...)

function Base.iterate(s::CellTypeSet)
    for i in UInt8(0):UInt8(N_CELL_TYPES)
        if i in s
            return (i, i)
        end
    end
    return nothing
end
function Base.iterate(s::CellTypeSet, prev_value::UInt8)
    for i in (prev_value+UInt8(1)):UInt8(N_CELL_TYPES)
        if i in s
            return (i, i)
        end
    end
    return nothing
end


#####################################
##   Cell lookups by type

const CELL_CODE_BY_CHAR::Dict{Char, UInt8} = Dict(
    t.char => t.code for t in CELL_TYPES
)

"Lookup by Symbol is optimal when parsing/interpreting pragmas"
const CELL_CODE_BY_SYMBOL::Dict{Symbol, UInt8} = Dict(Symbol(k)=>v for (k,v) in CELL_CODE_BY_CHAR)

"Represents 'null' or an unset cell"
const CELL_CODE_INVALID = convert(UInt8, 255)
"Represents 'wildcard' or 'unchanged'"
const CELL_CHAR_INVALID = '_'

CellTypeSet(chars::Char...) = CellTypeSet(getindex.(Ref(CELL_CODE_BY_CHAR), chars))
function CellTypeSet(chars::AbstractVector{Char})
    s = CellTypeSet()
    for c in chars
        s = push(s, CELL_CODE_BY_CHAR[c])
    end
    return s
end
function CellTypeSet(chars::String)
    s = CellTypeSet()
    for c in chars
        s = push(s, CELL_CODE_BY_CHAR[c])
    end
    return s
end
CellTypeSet(s::Symbol) = CellTypeSet(string(s))


#####################################
##   Cell grid

"An abstract grid of cell values (represented by their code)"
const CellGrid{N} = AbstractArray{UInt8, N}
"A contiguous heap-allocated grid of cell values (represented by their code)"
const CellGridConcrete{N} = Array{UInt8, N}

"Reference to a N-dimensional grid cell"
const CellIdx{N} = Vec{N, Int32}


"
An orthogonal direction along some N-dimensional grid.
The set of all possible directions is well-ordered
    (1/-1, 1/+1, 2/-1, 2/+1, etc),
    so you can refer to them by index as well.
"
struct GridDir
    # 1 - N
    axis::Int32
    # -1 or +1
    sign::Int32

    GridDir(axis, sign) = if (axis < 1) || !in(sign, (-1, +1))
        error("Axis must be >0 and sign must be -1 or +1! Got ", axis, " and ", sign)
    else
        new(convert(Int32, axis), convert(Int32, sign))
    end
    GridDir(axis, sign, skip_check) = if skip_check
        new(convert(Int32, axis), convert(Int32, sign))
    else
        GridDir(axis, sign)
    end
end

"Gets the index of the given grid direction"
grid_dir_index(d::GridDir)::Int32 = (d.axis * Int32(2)) + ((d.sign + 1) ÷ 2)
"Gets the grid direction at the given index"
grid_dir_index(_i) = let i = convert(Int32, _i)
    GridDir((i + 1) ÷ 2, (2 * ((i - 1) % 2)) - 1)
end
"Gets the number of grid directions for the given number of dimensions"
@inline grid_dir_count(N) = N * convert(typeof(N), 2)

function grid_dir_pos_along(dir::GridDir, pos::CellIdx, offset_along_dir::Integer)::typeof(pos)
    out_pos = pos
    @set! out_pos[dir.axis] += dir.sign * convert(Int32, offset_along_dir)
    return out_pos
end

# Optimize hashing by compressing the GridDir into its index.
Base.hash(d::GridDir, u::UInt) = Base.hash(grid_dir_index(d), u)

"A contiguous and oriented line of cells; see `CellRegion{N}` for the multidimensional version"
struct CellLine{N}
    start_cell::CellIdx{N}
    movement::GridDir
    length::Int32
end
function Base.print(io::IO, l::CellLine)
    p_start = l.start_cell
    p_end = p_start
    @set! p_end[l.movement.axis] += l.movement.sign * (l.length - 1)
    print(io, "<", p_start, " to ", p_end, ">")
end

"
A contiguous region of cells, the multidimensional equivalent of `CellLine{N}`.
It's just a wrapper around `BoxI{N}` but the goal is to be a clearly-typed sibling of `CellLine`.
"
struct CellRegion{N}
    b::BoxI{N}
end
CellRegion(r::Bplus.Math.VecRange) = CellRegion(Box(r))
Base.print(io::IO, r::CellRegion) = print(io, "R{", r.b, "}")

"
Walks through every cell in a small contiguous line/area.
Invokes your lambda with `local_idx, cell_idx`, where `local_idx` is either
  `Int32` for strips (`CellLine` input) or `VecI` for multidimensional (`CellRegion` input).

If you pass `true` then you can break the loop early by making your lambda return `true`,
  in which case this returns whether an early break happened.
"
function for_each_cell(toDo, area::Union{CellLine{N}, CellRegion{N}},
                       ::Val{Breakable}
                      )::(Breakable ? Bool : Nothing) where {N, Breakable}
    if area isa CellLine{N}
        for iz in zero(Int32):(area.length - one(Int32))
            cell = area.start_cell
            @set! cell[area.movement.axis] += (iz * area.movement.sign)
            result = toDo(iz + one(Int32), cell)
            if Breakable && convert(Bool, result)
                return true
            end
        end
    elseif area isa CellRegion{N}
        for iv in area.b
            result = toDo(iv - min_inclusive(area.b) + one(Int32), iv)
            if Breakable && convert(Bool, result)
                return true
            end
        end
    else
        error("Unhandled: ", typeof(area))
    end
    return Breakable ? false : nothing
end
@inline for_each_cell(toDo, line, breakable::Bool=false) = for_each_cell(toDo, line, Val(breakable))

"Gets the bounds around a `CellLine` or `CellRegion` (in the latter case it just returns the region's Box)"
function cell_line_aabb(area::Union{CellLine{N}, CellRegion{N}})::Box{N, Int32} where {N}
    if area isa CellRegion
        return area.b
    elseif area isa CellLine
        start_pos = area.start_cell
        end_pos = start_pos
        @set! end_pos[area.movement.axis] += area.movement.dir * (area.length - 1)
        (start_pos, end_pos) = Bplus.Math.minmax(start_pos, end_pos)
        return Box{N, Int32}(
            min=start_pos,
            max=end_pos
        )
    else
        error("Unhandled: ", typeof(area))
    end
end