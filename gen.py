import os

def getName(st):
    return st[0:st.find(".")]
def change(st):
    rt = ""
    for i in st:
        if st == '-':
            rt += '_'
        elif i.isupper(): 
            rt += i.lower()
        elif i.islower() or i == '_':
            rt += i
    return rt

codes = []
def write_into_markdown(root, dirs, files):
    path = os.path.join(root, dirs, files)
    dir_path = os.path.join('_code', dirs)
    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
    target = os.path.join("_code", dirs, getName(files) + ".md")
    # print(path, target)
    with open(path, "r") as fin:
        with open(target, "w") as fout:
            content = fin.read()
            fout.write('```cpp\n')
            try:
                fout.write(content)
            except e:
                fout.write(content.decode('gb2312').decode('utf-8'))
            fout.write('\n```\n')
            codes.append((files, content))
    return target
    
def walk_code(fp):
    root = "code"
    for dirs in os.listdir(root):
        path = os.path.join(root, dirs)
        if os.path.isdir(path):
            fp.write("- %s\n" % dirs)
            for files in os.listdir(path):
                file_path = os.path.join(path, files)
                if os.path.isfile(file_path):
                    fp.write("  - [%s](#%s)\n" % 
                            (
                                getName(files),
                                change(getName(files))
                            ))
                    write_into_markdown(root, dirs, files)

with open("codelist.md", "w") as fp:
    walk_code(fp)
    for (name, content) in codes:
        fp.write("\n\n")
        fp.write("# %s\n\n" % change(getName(name)))
        fp.write("```cpp\n")
        fp.write(content)
        fp.write("\n```\n")
