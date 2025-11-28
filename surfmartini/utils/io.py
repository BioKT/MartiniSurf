import os


def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)


def write_list(lst, handle, chunk=15):
    for i in range(0, len(lst), chunk):
        handle.write(" ".join(str(x) for x in lst[i : i + chunk]) + "\n")
