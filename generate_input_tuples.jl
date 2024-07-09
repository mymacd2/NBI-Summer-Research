
logαs = range(-26.0, -24.0, 11)
αs = 10 .^ logαs

# nubars = [true, false]

progs = ["small", "large", "bh"]

# cases = ["B", "C", "democratic"]

function generate_tuples()
    open("input_pairs.txt", "w") do file
        for sm in progs
            for α1 in αs
                for α2 in αs
                    println(file, "$bar $sm $α1 $α2")
                end
            end
        end
    end
end

generate_tuples()
