logαs = range(-27, -23, 21)
αs = 10 .^ logαs

daughters = ["1", "2", "alt"]

progs = ["small", "large", "bh"]

function generate_tuples()
    open("input_tuples_2nu.txt", "w") do file
        for daughter in daughters
            for prog in progs
                for α in αs
                    println(file, "$daughter $prog $α")
                end
            end
        end
    end
end

generate_tuples()