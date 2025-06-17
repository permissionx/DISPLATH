module Recorder     
export @record   
const FLUSH_BYTES = 8_192           

const _fh  = Dict{Symbol, IO}()     
const _buf = Dict{Symbol, IOBuffer}()

function _ensure(sym::Symbol)
    haskey(_fh, sym) && return
    io = open(string(sym), "w")    
    _fh[sym]  = io
    _buf[sym] = IOBuffer()
end

function _flush!(sym::Symbol)
    buf = _buf[sym]; io = _fh[sym]
    write(io, take!(buf)); flush(io)
end

macro record(file, value)
    quote
        local _sym = $(QuoteNode(file))    
        Recorder._ensure(_sym)
        local buf = Recorder._buf[_sym]
        print(buf, $(esc(value)), '\n')
        if Recorder.position(buf) ≥ Recorder.FLUSH_BYTES
            Recorder._flush!(_sym)
        end
    end
end

atexit() do
    for s in keys(_fh)
        _flush!(s); close(_fh[s])
    end
end

end

using .Recorder