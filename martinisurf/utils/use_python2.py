def find_python2():
    import shutil

    # Try common python2 names
    candidates = ["python2.7", "python2", "/usr/bin/python2.7", "/usr/bin/python2"]

    for exe in candidates:
        if shutil.which(exe):
            return exe

    raise RuntimeError(
        "❌ DNA mode requires python2.7, but no python2 interpreter was found.\n"
        "Install it with:\n"
        "    conda install -n martinisurf-dna python=2.7\n"
    )
