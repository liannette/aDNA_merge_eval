with open("Vi33.19_max250.txt", "w") as outfile:
    for line in outfile:
        num = int(line.rstrip())
        if num > 250:
            num = 250
        print(num, file=outfile)
