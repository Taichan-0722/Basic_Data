number = ["", "", ""]

with open("coordinate_xtd.txt") as file:
    for i, string in enumerate(file):
        for char in string:
            if char.isnumeric() or char == "." or char == "-":
                number[i] += char
        if ".." in number[i]:
            number[i] = float(number[i].strip(".."))
        else:
            number[i] = float(number[i].strip("."))

print(number[0], number[1], number[2])