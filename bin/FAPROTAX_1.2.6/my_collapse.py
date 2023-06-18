with open("./FAPROTAX.txt", "r") as f:
    data = f.readlines()

sep = "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"


def get_tax(line):
    return line.split("\t")[0].replace("*", " ").strip()


def get_group(data, index):
    while data[index].startswith("#"):
        index -= 1
    return data[index].split("\t")[0]


members: dict = {}
for index, line in enumerate(data):
    line = line.strip()
    if sep in line:
        group = get_group(data, index)
        members[group] = set()
        continue
    if line.startswith("*"):
        members[group].add(get_tax(line))
