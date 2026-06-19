function read_dat(filename::String; header_lines::Int=1)
  lines = readlines(filename)
  len = length(lines) - header_lines

  points = Array{Float64, 2}(undef, len, 2)

  for i in 1:len
    line = lines[header_lines + i]
    coords = split(line)
    points[i, 1] = parse(Float64, coords[1])
    points[i, 2] = parse(Float64, coords[2])
  end

  return points
end
