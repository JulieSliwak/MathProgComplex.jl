"""
max_err, errors = compare_dat(file1, file2, epsilon=1e-10)

Compute max_err, the maximum difference between common lines of both files, and
a dict associating an error to the line(s) where it occured.
"""
function compare_dat(file1::String, file2::String; epsilon = 1e-10, display_level=0)
  temp1, temp2 = "sorted_dat1.temp", "sorted_dat2.temp"

  sort_dat(file1, temp1)
  sort_dat(file2, temp2)

  data1, data2 = readdlm(temp1), readdlm(temp2)

  rm(temp1)
  rm(temp2)

  lines1 = SortedSet([data1[i, 1]*data1[i, 2]*data1[i, 3]*data1[i, 4] for i=1:size(data1, 1)])
  lines2 = SortedSet([data2[i, 1]*data2[i, 2]*data2[i, 3]*data2[i, 4] for i=1:size(data2, 1)])

  if size(data1, 1) != size(data2, 1)
    warn("dat files have different number of lines  ($(size(data1, 1)) and $(size(data2, 1)))")
  end

  if lines1 != lines2
    warn("dat files have different key entries\nfile1 \\ file2: $(length(setdiff(lines1, lines2)))\nfile2 \\ file1: $(length(setdiff(lines2, lines1)))")
  end

  if display_level>0
    println("nfile1 \\ file2: $(sort(collect(setdiff(lines1, lines2)))) \nfile2 \\ file1: $(sort(collect(setdiff(lines2, lines1))))")
  end

  errors = SortedDict{Float64, Any}()
  for i=1:min(size(data1, 1), size(data2, 1))
    if data1[i, 1:4] == data2[i, 1:4]
      cur_error = max(abs(data1[i, 5]-data2[i, 5]), abs(data1[i, 6]-data2[i, 6]))
      if cur_error > epsilon
        if !haskey(errors, cur_error)
          errors[cur_error] = SortedSet()
        end
        push!(errors[cur_error], (i, data1[i, 1:4]))
      end
    else
      # warn("Compare_dat(): line $i - mismatched on line entry")
    end
  end

  max_coeff_err = 0
  if length(errors) > 0
    max_coeff_err = maximum(collect(keys(errors)))
  end
  return max_coeff_err, errors
end


"""
sort_dat(filename, outfile)

Writes an outfile dat file with lexicographically sorted lines from filename.
"""
function sort_dat(filename, outfile)
  lines = SortedSet{String}()

  open(filename, "r") do f
    while !eof(f)
      line = readline(f)
      push!(lines, line)
    end
  end

  if isfile(outfile)
    rm(outfile)
  end

  open(outfile, "a") do f
    for line in sort(collect(lines))
      println(f, line)
    end
  end
end
