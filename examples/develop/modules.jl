function RandomInSquare(a::Float64, b::Float64)
    x = rand() * a
    y = rand() * b 
    return [x,y,0.0]
end