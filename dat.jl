# open("data.txt", "a") do f
    for _ in 1:9*(10^11)
        n1 = rand()
        n2 = rand()
        n3 = rand()
        write(f, "$n1 $n2 $n3 \n")
    end
end